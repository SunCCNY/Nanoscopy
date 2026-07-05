"""
ViT-PF.py

Permutation-Free Vision Transformer (ViT-PF) for SMLM emitter localization.

Treats localization as a CLASSIFICATION problem — a probability mass
distribution over all Jy*Jx subpixels of the FOV.

Architecture (standard ViT classification head)
------------------------------------------------
Input:  [CLS, p1, p2, ..., p49]  50 tokens (CLS + one per pixel)
        + learnable positional embedding (50 × d_model)
        → Multi-head self-attention × n_layers (standard Transformer encoder)
        → CLS token output (d_model)
        → Linear (d_model → Jy*Jx)
        → Softmax → probability map over Jy*Jx subpixels

For 7×7 pixels, r=10: Jy=Jx=70, output = 4900 probabilities.

Training
--------
Target: 2D Gaussian soft distribution centered at each true subpixel
        (σ=1.0 subpixel), normalized to sum to 1 over all Jy*Jx subpixels.
        Smoother than hard 1/M spike — helps separate nearby emitters.
Loss:   KL divergence with Gaussian soft targets.
No Hungarian assignment — subpixel assignment is O(M) quantization.
Frames where two emitters share a subpixel are discarded (~0.1%).

Inference
---------
Given M, pick top-M subpixels from the 4900-element softmax output.
Estimated location = centre of picked subpixel.
The full 70×70 probability map is a super-resolution image for biologists.

Advantages over ViT with Hungarian loss
----------------------------------------
1. No Hungarian assignment → training is O(M), not O(M³)
2. M-independent output — one model handles any M
3. Unknown M handled by probability thresholding
4. Probability map = super-resolution image for biologists

Reference
---------
Dosovitskiy et al. "An image is worth 16x16 words." ICLR 2021.
Sun, Y. (2025). Permutation-free SMLM localization via subpixel
classification with Vision Transformer.
"""

import math
import torch
import torch.nn as nn
import torch.nn.functional as F

from .config import Camera, Gaussian2DPSF, Noise, EmitterData, SystemConfig
from .position import sample_emitters_batch
from .Gauss2D import gauss2d_frame_torch


# ======================================================================
# Batch generator
# ======================================================================

def generate_batch_pf(
    batch_size: int,
    M:          int,
    camera:     Camera,
    psf:        Gaussian2DPSF,
    noise:      Noise,
    Im0:        float,
    r:          int,
    margin_px:  int = 2,
) -> tuple:
    """Generate frames and soft subpixel target distributions.

    Parameters
    ----------
    batch_size : int
    M          : number of emitters per frame
    camera     : Camera
    psf        : Gaussian2DPSF
    noise      : Noise
    Im0        : emitter intensity (photons/s)
    r          : integer upsampling factor (Jx=Kx*r, Jy=Ky*r)
    margin_px  : edge pixels excluded when sampling positions

    Returns
    -------
    x      : (B, Kx*Ky)  float32 — normalised flattened frames
    target : (B, Jy*Jx)  float32 — Gaussian soft target summing to 1
    xy_nm  : (B, M, 2)   float32 — true emitter positions in nm
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    dtype  = torch.float32

    Kx, Ky = camera.Kx, camera.Ky
    Dx, Dy = camera.Dx, camera.Dy
    Jx, Jy = Kx * r, Ky * r
    dx, dy = Dx / r, Dy / r

    x_out      = torch.zeros(batch_size, Kx * Ky, dtype=dtype, device=device)
    target_out = torch.zeros(batch_size, Jy * Jx, dtype=dtype, device=device)
    xy_out     = torch.zeros(batch_size, M, 2,    dtype=dtype, device=device)

    n_done = 0

    while n_done < batch_size:
        remain = batch_size - n_done
        xy = sample_emitters_batch(remain, M, camera, margin_px=margin_px)

        # Quantize to subpixel indices
        ix = torch.clamp((xy[:, :, 0] / dx).long(), 0, Jx - 1)  # (remain, M)
        iy = torch.clamp((xy[:, :, 1] / dy).long(), 0, Jy - 1)  # (remain, M)
        idx = iy * Jx + ix                                        # (remain, M)

        # Discard frames with collisions (two emitters in same subpixel)
        idx_sorted, _ = idx.sort(dim=1)
        collision = (idx_sorted[:, 1:] == idx_sorted[:, :-1]).any(dim=1)
        valid     = ~collision

        n_valid = valid.sum().item()
        if n_valid == 0:
            continue

        xy_v   = xy[valid][:batch_size - n_done]    # (n_take, M, 2)
        idx_v  = idx[valid][:batch_size - n_done]   # (n_take, M)
        n_take = xy_v.shape[0]

        # Generate frames
        Im_tensor = torch.full((n_take, M), float(Im0), dtype=dtype, device=device)
        emitter   = EmitterData(xy=xy_v, Im=Im_tensor)
        system    = SystemConfig(emitter=emitter, psf=psf, camera=camera, noise=noise)
        with torch.no_grad():
            frames = gauss2d_frame_torch(system)   # (n_take, Ky, Kx)

        # Normalise frames
        x_flat = frames.view(n_take, -1)
        x_flat = x_flat / (x_flat.amax(dim=1, keepdim=True) + 1e-8)

        # Hard soft target: 1/M at each true subpixel, 0 elsewhere
        tgt = torch.zeros(n_take, Jy * Jx, dtype=dtype, device=device)
        b_idx = torch.arange(n_take, device=device).unsqueeze(1).expand(n_take, M)
        tgt.scatter_(1, idx_v, 1.0 / M)

        x_out[n_done:n_done + n_take]      = x_flat
        target_out[n_done:n_done + n_take] = tgt
        xy_out[n_done:n_done + n_take]     = xy_v
        n_done += n_take

    return x_out, target_out, xy_out


# ======================================================================
# ViT-PF Localizer
# ======================================================================

class ViTPFLocalizer(nn.Module):
    """Permutation-Free ViT for SMLM — standard CLS classification head.

    Standard ViT encoder with CLS token, followed by a single linear
    layer mapping to Jy*Jx subpixel logits. Softmax produces the
    probability mass distribution over all subpixels.

    Input:  (B, Kx*Ky) normalised flattened frame
    Output: (B, Jy*Jx) subpixel logits (softmax → probability map)

    Parameters
    ----------
    Kx, Ky   : frame dimensions in pixels
    r        : upsampling factor (subpixel size = Dx/r nm, default 10)
    d_model  : token embedding dimension (default 128)
    n_heads  : attention heads (default 4)
    n_layers : transformer encoder blocks (default 6)
    d_ff     : feed-forward hidden dim (default 512)
    dropout  : dropout probability (default 0.1)
    """

    def __init__(
        self,
        Kx:      int,
        Ky:      int,
        r:       int   = 10,
        d_model: int   = 256,
        n_heads: int   = 8,
        n_layers:int   = 6,
        d_ff:    int   = 512,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.Kx      = Kx
        self.Ky      = Ky
        self.r       = r
        self.Jx      = Kx * r
        self.Jy      = Ky * r
        self.d_model = d_model
        n_tokens     = Kx * Ky   # one token per pixel

        # Project each pixel scalar → d_model
        self.patch_embed = nn.Linear(1, d_model)

        # Learnable CLS token and positional embeddings
        self.cls_token = nn.Parameter(torch.zeros(1, 1, d_model))
        self.pos_embed = nn.Parameter(
            torch.zeros(1, n_tokens + 1, d_model)   # +1 for CLS
        )
        nn.init.trunc_normal_(self.cls_token, std=0.02)
        nn.init.trunc_normal_(self.pos_embed, std=0.02)

        self.drop = nn.Dropout(dropout)

        # Standard Transformer encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model         = d_model,
            nhead           = n_heads,
            dim_feedforward = d_ff,
            dropout         = dropout,
            activation      = 'gelu',
            batch_first     = True,
            norm_first      = False,
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=n_layers)

        # Joint head: ALL token representations → Jy*Jx subpixel logits
        # Input: (n_tokens+1) * d_model = 50 * d_model
        # Every subpixel logit is computed from the full frame representation
        n_flat = (Kx * Ky + 1) * d_model   # 50 * d_model
        self.head = nn.Sequential(
            nn.Linear(n_flat, d_ff),
            nn.GELU(),
            nn.Linear(d_ff, self.Jy * self.Jx),
        )

        self._init_weights()

    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.trunc_normal_(m.weight, std=0.02)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)
            elif isinstance(m, nn.LayerNorm):
                nn.init.ones_(m.weight)
                nn.init.zeros_(m.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        x : (B, Kx*Ky) normalised flattened frame
        returns: (B, Jy*Jx) subpixel logits (apply softmax for probabilities)
        """
        B = x.shape[0]

        # Embed pixels → tokens
        tokens = self.patch_embed(x.unsqueeze(-1))    # (B, n_tokens, d_model)

        # Prepend CLS token
        cls    = self.cls_token.expand(B, -1, -1)     # (B, 1, d_model)
        tokens = torch.cat([cls, tokens], dim=1)       # (B, n_tokens+1, d_model)

        # Add positional embedding
        tokens = self.drop(tokens + self.pos_embed)

        # Transformer encoder — all 50 tokens attend to each other
        tokens = self.encoder(tokens)                  # (B, n_tokens+1, d_model)

        # Flatten ALL 50 tokens (CLS + 49 pixel tokens) into one vector
        # (B, n_tokens+1, d_model) → (B, (n_tokens+1)*d_model)
        flat = tokens.reshape(B, -1)                   # (B, 50*d_model)

        # Joint linear head: all token information → 1225 subpixel logits
        return self.head(flat)                         # (B, Jy*Jx)

    def predict_xy(
        self,
        x:          torch.Tensor,
        M:          int,
        camera:     Camera,
        jitter:     bool  = False,
        nms_radius: float = 30.0,
    ) -> torch.Tensor:
        """Inference: sequential peak detection with Non-Maximum Suppression (NMS).

        For m = 1, ..., M:
          (i)  Pick peak subpixel s' = argmax p^(m)(s)
               Set estimated location mu_m = centre of s'
          (ii) Zero out all subpixels within nms_radius nm of mu_m:
               p^(m+1)(s) = 0             if ||c_s - mu_m|| <= nms_radius
               p^(m+1)(s) = p^(m)(s)      otherwise

        Parameters
        ----------
        x          : (B, Kx*Ky) normalised flattened frames
        M          : number of emitters
        camera     : Camera
        jitter     : if True, add U(-dx/2,dx/2) x U(-dy/2,dy/2) jitter
        nms_radius : suppression radius in nm (default 30 nm)

        Returns
        -------
        xy_est : (B, M, 2) estimated emitter positions in nm
        """
        self.eval()
        dx = camera.Dx / self.r
        dy = camera.Dy / self.r
        B  = x.shape[0]
        device = x.device

        with torch.no_grad():
            logits = self.forward(x)               # (B, Jy*Jx)
            probs  = torch.softmax(logits, dim=1)  # (B, Jy*Jx)

        # Subpixel centre coordinates in nm
        gy = torch.arange(self.Jy, device=device, dtype=torch.float32)
        gx = torch.arange(self.Jx, device=device, dtype=torch.float32)
        cy_flat = (gy.unsqueeze(1).expand(self.Jy, self.Jx).reshape(-1) + 0.5) * dy
        cx_flat = (gx.unsqueeze(0).expand(self.Jy, self.Jx).reshape(-1) + 0.5) * dx
        cx_flat = cx_flat.unsqueeze(0)   # (1, S)
        cy_flat = cy_flat.unsqueeze(0)   # (1, S)

        xy_est   = torch.zeros(B, M, 2, device=device)
        residual = probs.clone()          # (B, S)

        for m in range(M):
            # (i) Pick peak subpixel
            peak_idx = residual.argmax(dim=1)              # (B,)
            mu_x = cx_flat[0, peak_idx]                    # (B,)
            mu_y = cy_flat[0, peak_idx]                    # (B,)

            xy_est[:, m, 0] = mu_x
            xy_est[:, m, 1] = mu_y

            # (ii) NMS: zero out all subpixels within nms_radius of peak
            dist2    = (cx_flat - mu_x.unsqueeze(1))**2 + \
                       (cy_flat - mu_y.unsqueeze(1))**2    # (B, S)
            suppress = dist2 <= nms_radius**2              # (B, S) boolean mask
            residual = residual * (~suppress).float()      # zero suppressed region

        # Optional uniform jitter
        if jitter:
            xy_est[..., 0] += (torch.rand_like(xy_est[..., 0]) - 0.5) * dx
            xy_est[..., 1] += (torch.rand_like(xy_est[..., 1]) - 0.5) * dy

        return xy_est   # (B, M, 2)


    def predict_xy_from_probmap(
        self,
        prob_map:   torch.Tensor,
        M:          int,
        camera:     Camera,
        jitter:     bool  = False,
        nms_radius: float = 30.0,
    ) -> torch.Tensor:
        """Sequential peak detection directly on a given probability map.

        Same algorithm as predict_xy but operates on a pre-computed
        probability map instead of raw frames. Used for averaged maps.

        Parameters
        ----------
        prob_map : (Jy, Jx) or (B, Jy, Jx) probability map
        M        : number of emitters
        camera   : Camera
        jitter   : add uniform random offset within subpixel

        Returns
        -------
        xy_est : (M, 2) or (B, M, 2) estimated positions in nm
        """
        dx = camera.Dx / self.r
        dy = camera.Dy / self.r
        device = prob_map.device

        single = prob_map.ndim == 2
        if single:
            prob_map = prob_map.unsqueeze(0)   # (1, Jy, Jx)
        B = prob_map.shape[0]

        # Flatten to (B, Jy*Jx)
        probs = prob_map.view(B, -1)

        # Subpixel centre coordinates
        gy = torch.arange(self.Jy, device=device, dtype=torch.float32)
        gx = torch.arange(self.Jx, device=device, dtype=torch.float32)
        cy_flat = (gy.unsqueeze(1).expand(self.Jy, self.Jx).reshape(-1) + 0.5) * dy
        cx_flat = (gx.unsqueeze(0).expand(self.Jy, self.Jx).reshape(-1) + 0.5) * dx
        cx_flat = cx_flat.unsqueeze(0)   # (1, S)
        cy_flat = cy_flat.unsqueeze(0)   # (1, S)

        xy_est   = torch.zeros(B, M, 2, device=device)
        residual = prob_map.reshape(B, -1).clone()   # (B, S) flattened

        for m in range(M):
            # (i) Pick peak subpixel
            peak_idx = residual.argmax(dim=1)
            mu_x = cx_flat[0, peak_idx]
            mu_y = cy_flat[0, peak_idx]

            xy_est[:, m, 0] = mu_x
            xy_est[:, m, 1] = mu_y

            # (ii) NMS: zero out all subpixels within nms_radius of peak
            dist2    = (cx_flat - mu_x.unsqueeze(1))**2 + \
                       (cy_flat - mu_y.unsqueeze(1))**2
            suppress = dist2 <= nms_radius**2
            residual = residual * (~suppress).float()

        if jitter:
            xy_est[..., 0] += (torch.rand_like(xy_est[..., 0]) - 0.5) * dx
            xy_est[..., 1] += (torch.rand_like(xy_est[..., 1]) - 0.5) * dy

        return xy_est.squeeze(0) if single else xy_est   # (M,2) or (B,M,2)

    def predict_prob_map(
        self,
        x: torch.Tensor,
    ) -> torch.Tensor:
        """Return (Jy, Jx) probability map — super-resolution image."""
        single = x.ndim == 1
        if single:
            x = x.unsqueeze(0)
        B = x.shape[0]

        self.eval()
        with torch.no_grad():
            logits   = self.forward(x)
            probs    = torch.softmax(logits, dim=1)    # (B, Jy*Jx)
            prob_map = probs.view(B, self.Jy, self.Jx)

        return prob_map.squeeze(0) if single else prob_map
