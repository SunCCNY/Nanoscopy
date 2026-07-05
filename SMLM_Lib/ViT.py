"""
ViT.py

Vision Transformer (ViT) for SMLM emitter position estimation.

Each 7×7 frame is treated as a sequence of 49 tokens (one per pixel).
A learnable CLS token aggregates global information; the transformer
encoder applies multi-head self-attention over all 50 tokens; the CLS
token output is passed through an MLP regression head to produce M
emitter positions in nm.

Uses PyTorch built-in nn.TransformerEncoderLayer with batch_first=True
for optimised fused kernels, and torch.compile() for kernel fusion.

Loss: hungarian_loss from train_utils.py (permutation-invariant, 1-to-1 matching).
Batch generator: generate_batch from train_utils.py (on-the-fly, no disk I/O).

Reference
---------
Dosovitskiy, A. et al. "An image is worth 16x16 words: Transformers for
image recognition at scale." ICLR 2021.

Adapted for single-molecule localization by Yi Sun (2025).
"""

import math
import torch
import torch.nn as nn
import torch.nn.functional as F

from .config import Camera, Gaussian2DPSF, Noise, EmitterData, SystemConfig
from .position import sample_emitters_batch
from .Gauss2D import gauss2d_frame_torch

# Re-export shared utilities so ViT-Training.py can import from one place
from .train_utils import generate_batch, hungarian_loss


# ======================================================================
# ViT Localizer
# ======================================================================

class ViTLocalizer(nn.Module):
    """Vision Transformer localizer for SMLM.

    Treats each pixel of a (Ky, Kx) frame as one token.  A learnable
    CLS token aggregates the sequence; its final representation is
    passed through an MLP head to regress M emitter positions in nm.

    Uses nn.TransformerEncoderLayer (batch_first=True) which employs
    PyTorch's optimised fused attention kernels and is compatible with
    torch.compile() for additional kernel-fusion speedup.

    Input:  (B, Kx*Ky) normalised flattened frame
    Output: (B, M, 2)  emitter positions in nm

    Parameters
    ----------
    Kx, Ky    : frame dimensions in pixels
    M         : number of emitters per frame
    d_model   : token embedding dimension (default 128)
    n_heads   : number of attention heads (default 4)
    n_layers  : number of transformer encoder blocks (default 6)
    d_ff      : feed-forward hidden dimension (default 512)
    dropout   : dropout probability (default 0.1)

    Notes
    -----
    For a 7×7 frame: 49 pixel tokens + 1 CLS token = 50 tokens.
    With d_model=128, n_heads=4, n_layers=6 this model has ~1.25M
    parameters — 6× CNN-KAN — providing sufficient capacity to match
    EM/UGIA-F accuracy at high emitter density.
    The smaller variant (d_model=32, n_heads=2, n_layers=2, ~32k params)
    is too small to converge for M≥2.
    """

    def __init__(
        self,
        Kx:       int,
        Ky:       int,
        M:        int,
        d_model:  int   = 128,
        n_heads:  int   = 4,
        n_layers: int   = 6,
        d_ff:     int   = 512,
        dropout:  float = 0.1,
    ):
        super().__init__()
        self.Kx      = Kx
        self.Ky      = Ky
        self.M       = M
        self.d_model = d_model
        n_tokens     = Kx * Ky   # one token per pixel

        # Project each pixel (scalar) → d_model
        self.patch_embed = nn.Linear(1, d_model)

        # Learnable CLS token and positional embeddings
        self.cls_token = nn.Parameter(torch.zeros(1, 1, d_model))
        self.pos_embed = nn.Parameter(
            torch.zeros(1, n_tokens + 1, d_model)   # +1 for CLS
        )
        nn.init.trunc_normal_(self.cls_token, std=0.02)
        nn.init.trunc_normal_(self.pos_embed, std=0.02)

        self.drop = nn.Dropout(dropout)

        # Transformer encoder — uses PyTorch built-in for fused kernels.
        # norm_first=False (default) enables the fast path optimized
        # implementation in nn.TransformerEncoderLayer.
        encoder_layer = nn.TransformerEncoderLayer(
            d_model         = d_model,
            nhead           = n_heads,
            dim_feedforward = d_ff,
            dropout         = dropout,
            activation      = 'gelu',
            batch_first     = True,   # (B, T, d_model) convention
            norm_first      = False,  # post-norm: enables PyTorch fast path
        )
        self.encoder = nn.TransformerEncoder(
            encoder_layer,
            num_layers = n_layers,
        )

        # MLP regression head: CLS token → M * 2 coordinates
        self.head = nn.Sequential(
            nn.Linear(d_model, d_ff),
            nn.GELU(),
            nn.Linear(d_ff, M * 2),
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
        returns: (B, M, 2) emitter positions in nm
        """
        B = x.shape[0]

        # Embed each pixel as a token: (B, n_tokens, 1) → (B, n_tokens, d_model)
        tokens = self.patch_embed(x.unsqueeze(-1))    # (B, n_tokens, d_model)

        # Prepend CLS token
        cls    = self.cls_token.expand(B, -1, -1)     # (B, 1, d_model)
        tokens = torch.cat([cls, tokens], dim=1)       # (B, n_tokens+1, d_model)

        # Add positional embedding and dropout
        tokens = self.drop(tokens + self.pos_embed)

        # Transformer encoder (fused PyTorch kernels)
        tokens = self.encoder(tokens)                  # (B, n_tokens+1, d_model)

        # CLS token → regression head
        cls_out = tokens[:, 0]                         # (B, d_model)
        return self.head(cls_out).view(B, self.M, 2)   # (B, M, 2)
