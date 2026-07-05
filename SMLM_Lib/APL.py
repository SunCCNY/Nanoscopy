"""
APL.py

Adaptive Piecewise-Linear (APL) localizer for SMLM emitter position
estimation.  This is the model previously labelled "KAN"/"CNN-KAN".

Each APL layer = linear mixing (W x + b) followed by a learnable
piecewise-linear activation (triangular-hat / first-order B-spline
basis on a learnable knot grid).  This is the adaptive-activation
design of Agostinelli et al. (2015) / Scardapane et al., NOT the
edge-function KAN of Liu et al. (2024) — see KAN.py for that.

Shared utilities generate_batch()/hungarian_loss() now live in
train_utils.py.

Yi Sun
"""

import torch
import torch.nn as nn
import torch.nn.functional as F

# ======================================================================
# APL Layer
# ======================================================================

class APLLayer(nn.Module):
    """Single APL layer with learnable B-spline-like basis functions.

    Parameters
    ----------
    in_dim  : input dimension
    out_dim : output dimension
    n_knots : number of basis knots (default 16)
    x_min   : left boundary of knot grid
    x_max   : right boundary of knot grid
    """

    def __init__(
        self,
        in_dim:  int,
        out_dim: int,
        n_knots: int   = 16,
        x_min:   float = -3.0,
        x_max:   float =  3.0,
    ):
        super().__init__()
        self.weight = nn.Parameter(torch.randn(out_dim, in_dim) * 0.1)
        self.bias   = nn.Parameter(torch.zeros(out_dim))
        self.knots  = nn.Parameter(torch.linspace(x_min, x_max, n_knots))
        self.coeffs = nn.Parameter(torch.randn(out_dim, n_knots) * 0.1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        z          = F.linear(x, self.weight, self.bias)
        z_expanded = z.unsqueeze(-1)
        knots      = self.knots.view(1, 1, -1)
        dist       = torch.abs(z_expanded - knots)

        with torch.no_grad():
            sorted_knots, _ = torch.sort(self.knots)
            diffs = sorted_knots[1:] - sorted_knots[:-1]
            h = max(diffs.mean().item(), 1e-3)

        basis = torch.clamp(1.0 - dist / h, min=0.0)
        return torch.sum(basis * self.coeffs.unsqueeze(0), dim=-1)


# ======================================================================
# APL Localizer
# ======================================================================

class APLLocalizer(nn.Module):
    """Pure APL localizer (no CNN frontend).

    Input:  (B, Kx*Ky) normalised flattened frame
    Output: (B, M, 2)  emitter positions in nm
    """

    def __init__(
        self,
        Kx:      int,
        Ky:      int,
        M:       int,
        hidden1: int = 256,
        hidden2: int = 128,
        n_knots: int = 16,
    ):
        super().__init__()
        self.M    = M
        in_dim    = Kx * Ky
        self.apl1 = APLLayer(in_dim,  hidden1, n_knots=n_knots)
        self.apl2 = APLLayer(hidden1, hidden2, n_knots=n_knots)
        self.out  = nn.Linear(hidden2, M * 2)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        B = x.shape[0]
        h = F.relu(self.apl1(x))
        h = F.relu(self.apl2(h))
        return self.out(h).view(B, self.M, 2)   # (B, M, 2)


# ======================================================================
# CNN-APL Localizer
# ======================================================================

class CNNAPLLocalizer(nn.Module):
    """CNN front-end + APL regression head for SMLM localization.

    Input:  (B, Kx*Ky) normalised flattened frame
    Output: (B, M, 2)  emitter positions in nm

    Default hyperparameters (n_channels=16, n_knots=16) are optimal
    for 7×7 frames — empirically verified.
    """

    def __init__(
        self,
        Kx:         int,
        Ky:         int,
        M:          int,
        n_channels: int = 16,
        hidden_dim: int = 256,
        n_knots:    int = 16,
    ):
        super().__init__()
        self.Kx = Kx
        self.Ky = Ky
        self.M  = M

        self.cnn = nn.Sequential(
            nn.Conv2d(1, n_channels, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
            nn.Conv2d(n_channels, n_channels, kernel_size=3, padding=1),
            nn.ReLU(inplace=True),
        )

        feat_dim  = n_channels * Kx * Ky
        self.apl1 = APLLayer(feat_dim,   hidden_dim, n_knots=n_knots)
        self.apl2 = APLLayer(hidden_dim, hidden_dim, n_knots=n_knots)
        self.out  = nn.Linear(hidden_dim, M * 2)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        B     = x.shape[0]
        x_img = x.view(B, 1, self.Ky, self.Kx)
        f     = self.cnn(x_img).view(B, -1)
        h     = self.apl1(f)
        h     = self.apl2(h)
        return self.out(h).view(B, self.M, 2)   # (B, M, 2)
