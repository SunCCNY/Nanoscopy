"""
KAN.py  (TRUE Kolmogorov-Arnold Network, Liu et al. 2024)

Faithful implementation of the KAN layer from
    Liu et al., "KAN: Kolmogorov-Arnold Networks", arXiv:2404.19756 (2024),
using the standard efficient B-spline formulation:

    phi(x) = w_base * SiLU(x)  +  w_spline * sum_j c_j * B_j(x)

  * NO linear weight matrix in the MLP sense.
  * Learnable univariate function on every (input, output) EDGE.
  * Each edge function = SiLU residual base + a B-spline of order k on a
    uniform grid of G intervals.
  * Node value = sum over input edges (summation at nodes).

This is the real KAN, in contrast to APL.py (linear mixing + learnable
piecewise-linear activation), which is what the earlier "KAN"/"CNN-KAN"
models actually were.

Drop-in companions:
    generate_batch(), hungarian_loss()  -> train_utils.py (shared)

Yi Sun
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


# ======================================================================
# True KAN layer (B-spline edges + SiLU residual)
# ======================================================================

class KANLayer(nn.Module):
    """A single Kolmogorov-Arnold layer with B-spline edge functions.

    Parameters
    ----------
    in_features  : input dimension  I
    out_features : output dimension O
    grid_size    : number of spline intervals G          (Liu default 5)
    spline_order : B-spline order k                       (Liu default 3, cubic)
    grid_range   : (lo, hi) range over which the grid is laid out
    scale_noise  : std of the noise used to initialise spline coefficients
    base_activation : residual base b(x)                  (Liu default SiLU)

    The edge function on edge (o, i) is
        phi_{o,i}(x) = base_weight[o,i] * b(x)
                     + spline_weight[o,i,:] . B(x)
    and the node output is  y_o = sum_i phi_{o,i}(x_i).
    """

    def __init__(
        self,
        in_features:  int,
        out_features: int,
        grid_size:    int   = 5,
        spline_order: int   = 3,
        grid_range:   tuple = (-1.0, 1.0),
        scale_noise:  float = 0.1,
        base_activation=nn.SiLU,
    ):
        super().__init__()
        self.in_features  = in_features
        self.out_features = out_features
        self.grid_size    = grid_size
        self.spline_order = spline_order

        # ----- fixed uniform grid, extended by `spline_order` on each side ----
        # shape: (in_features, grid_size + 2*spline_order + 1)
        h = (grid_range[1] - grid_range[0]) / grid_size
        grid = (
            torch.arange(-spline_order, grid_size + spline_order + 1) * h
            + grid_range[0]
        )
        grid = grid.expand(in_features, -1).contiguous()
        self.register_buffer("grid", grid)   # not learned (vanilla KAN, fixed grid)

        # ----- learnable parameters --------------------------------------
        # residual base scale, one per edge
        self.base_weight = nn.Parameter(
            torch.empty(out_features, in_features)
        )
        # spline coefficients: (out, in, n_basis) with n_basis = G + k
        self.spline_weight = nn.Parameter(
            torch.empty(out_features, in_features, grid_size + spline_order)
        )
        self.base_activation = base_activation()

        self.reset_parameters(scale_noise)

    # ------------------------------------------------------------------
    def reset_parameters(self, scale_noise: float = 0.1):
        nn.init.kaiming_uniform_(self.base_weight, a=5 ** 0.5)
        with torch.no_grad():
            # initialise spline coeffs from small noise sampled at grid points
            noise = (
                (torch.rand(self.grid_size + 1, self.in_features, self.out_features)
                 - 0.5) * scale_noise / self.grid_size
            )
            self.spline_weight.data.copy_(
                self._curve2coeff(
                    self.grid.T[self.spline_order:-self.spline_order],  # (G+1, I)
                    noise,
                )
            )

    # ------------------------------------------------------------------
    def b_splines(self, x: torch.Tensor) -> torch.Tensor:
        """Evaluate B-spline bases on the inputs.

        x : (batch, in_features)
        returns (batch, in_features, grid_size + spline_order)
        via the Cox-de Boor recursion.
        """
        grid = self.grid                       # (in, G + 2k + 1)
        x = x.unsqueeze(-1)                    # (batch, in, 1)
        bases = ((x >= grid[:, :-1]) & (x < grid[:, 1:])).to(x.dtype)
        for k in range(1, self.spline_order + 1):
            left = (x - grid[:, : -(k + 1)]) / (grid[:, k:-1] - grid[:, : -(k + 1)])
            right = (grid[:, k + 1:] - x) / (grid[:, k + 1:] - grid[:, 1:-k])
            bases = left * bases[:, :, :-1] + right * bases[:, :, 1:]
        return bases.contiguous()              # (batch, in, G + k)

    # ------------------------------------------------------------------
    def _curve2coeff(self, x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
        """Least-squares fit of spline coefficients to (x, y) samples.

        x : (n_samples, in_features)
        y : (n_samples, in_features, out_features)
        returns (out_features, in_features, grid_size + spline_order)
        """
        A = self.b_splines(x).transpose(0, 1)          # (in, n_samples, G+k)
        B = y.transpose(0, 1)                          # (in, n_samples, out)
        sol = torch.linalg.lstsq(A, B).solution        # (in, G+k, out)
        return sol.permute(2, 0, 1).contiguous()       # (out, in, G+k)

    # ------------------------------------------------------------------
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # residual base: SiLU(x) mixed by base_weight  -> (batch, out)
        base = F.linear(self.base_activation(x), self.base_weight)
        # spline part: bases (batch, in, G+k) contracted with spline_weight
        bases = self.b_splines(x).view(x.size(0), -1)              # (batch, in*(G+k))
        spline = F.linear(bases, self.spline_weight.view(self.out_features, -1))
        return base + spline

    # ------------------------------------------------------------------
    def regularization_loss(self, reg_activation=1.0, reg_entropy=1.0):
        """Optional L1 + entropy penalty on spline edges (Liu et al. sec. 2.5).

        Only needed if you want the sparsity / interpretability story;
        not required for pure localisation accuracy.
        """
        l1 = self.spline_weight.abs().mean(-1)        # (out, in)
        total = l1.sum()
        p = l1 / (total + 1e-12)
        entropy = -(p * (p + 1e-12).log()).sum()
        return reg_activation * total + reg_entropy * entropy


# ======================================================================
# Pure KAN localizer
# ======================================================================

class KANLocalizer(nn.Module):
    """Pure KAN localizer (Liu-faithful), no CNN front-end.

    Input:  (B, Kx*Ky) normalised flattened frame
    Output: (B, M, 2)  emitter positions in nm

    KANs are more expressive per edge than MLP neurons, so the hidden
    widths here are deliberately smaller than the APL baseline.  Treat
    (hidden, grid_size, spline_order) as the quantities to sweep.
    """

    def __init__(
        self,
        Kx:           int,
        Ky:           int,
        M:            int,
        hidden:       tuple = (64, 32),
        grid_size:    int   = 5,
        spline_order: int   = 3,
        layernorm:    bool  = True,
    ):
        super().__init__()
        self.M = M
        dims = [Kx * Ky, *hidden, M * 2]

        layers = []
        for i in range(len(dims) - 1):
            if layernorm and i > 0:
                layers.append(nn.LayerNorm(dims[i]))   # keep activations in grid range
            layers.append(
                KANLayer(dims[i], dims[i + 1],
                         grid_size=grid_size, spline_order=spline_order)
            )
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        B = x.shape[0]
        return self.net(x).view(B, self.M, 2)


# ======================================================================
# CNN + KAN localizer  (convolution front-end + true KAN head)
# ======================================================================

class CNNKANLocalizer(nn.Module):
    """Convolutional front-end + true-KAN regression head.

    Same conv stem idea as the APL version, feeding a Liu-faithful KAN
    head.  An optional 2x2 max-pool after the stem lets you *shrink* the
    head's input (and so the overall model) instead of expanding it.
    """

    def __init__(
        self,
        Kx:           int,
        Ky:           int,
        M:            int,
        n_channels:   int   = 16,
        hidden:       tuple = (64, 32),
        grid_size:    int   = 5,
        spline_order: int   = 3,
        pool:         bool  = False,
        layernorm:    bool  = True,
    ):
        super().__init__()
        self.Kx, self.Ky, self.M = Kx, Ky, M

        stem = [
            nn.Conv2d(1, n_channels, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv2d(n_channels, n_channels, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
        ]
        if pool:
            stem.append(nn.MaxPool2d(2))
        self.cnn = nn.Sequential(*stem)

        # infer flattened feature dim
        with torch.no_grad():
            dummy = torch.zeros(1, 1, Ky, Kx)
            feat_dim = self.cnn(dummy).view(1, -1).shape[1]

        dims = [feat_dim, *hidden, M * 2]
        layers = []
        for i in range(len(dims) - 1):
            if layernorm and i > 0:
                layers.append(nn.LayerNorm(dims[i]))
            layers.append(
                KANLayer(dims[i], dims[i + 1],
                         grid_size=grid_size, spline_order=spline_order)
            )
        self.head = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        B = x.shape[0]
        x_img = x.view(B, 1, self.Ky, self.Kx)
        f = self.cnn(x_img).view(B, -1)
        return self.head(f).view(B, self.M, 2)
