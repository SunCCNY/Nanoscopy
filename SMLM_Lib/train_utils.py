"""
train_utils.py

Shared training utilities for all learned SMLM localizers
(APL, CNN-APL, KAN, CNN-KAN, ViT):

  generate_batch()   — on-the-fly batch generator for training
  hungarian_loss()   — permutation-invariant loss via optimal 1-to-1 matching.
                       Solver selected automatically (fastest available):
                         1. torch_linear_assignment  — GPU-native, any M
                         2. Brute-force permutation  — pure PyTorch GPU, M <= 7
                         3. scipy                    — CPU fallback, any M

Moved here from KAN.py so the loss/generator no longer live inside a
model-specific module.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from itertools import permutations as _permutations

# Solver priority 1: GPU-native torch_linear_assignment.
# Test that the CUDA backend actually works (not just CPU fallback).
_SOLVER = 'scipy'   # default — overridden below if GPU works
_gpu_lap = None
try:
    from torch_linear_assignment import batch_linear_assignment as _gpu_lap_candidate
    # Verify the GPU backend is truly available by running a tiny test
    if torch.cuda.is_available():
        import warnings
        _test_C = torch.zeros(1, 2, 2, device='cuda')
        with warnings.catch_warnings():
            warnings.simplefilter("error")   # treat warning as error
            try:
                _gpu_lap_candidate(_test_C)
                _gpu_lap = _gpu_lap_candidate
                _SOLVER  = 'gpu'
            except Exception:
                pass   # CUDA backend not available — stay with scipy
    del _gpu_lap_candidate
except ImportError:
    pass

# Solver priority 3: scipy CPU fallback
from scipy.optimize import linear_sum_assignment as _scipy_lap

# Pre-compute all M! permutation tensors for M = 1..7 (stored on CPU,
# moved to device on first use).  M=7 has 5040 permutations — still fast.
_PERM_CACHE: dict = {}

def _get_perms(M: int, device: torch.device) -> torch.Tensor:
    """Return (M!, M) int64 tensor of all permutations of range(M)."""
    if M not in _PERM_CACHE:
        _PERM_CACHE[M] = torch.tensor(
            list(_permutations(range(M))), dtype=torch.long
        )
    return _PERM_CACHE[M].to(device)

from .config import Camera, Gaussian2DPSF, Noise, EmitterData, SystemConfig
from .position import sample_emitters_batch
from .Gauss2D import gauss2d_frame_torch


# ======================================================================
# Batch generator
# ======================================================================

def generate_batch(
    batch_size: int,
    M:          int,
    camera:     Camera,
    psf:        Gaussian2DPSF,
    noise:      Noise,
    Im0:        float = 300000.0,
    margin_px:  int   = 2,
) -> tuple:
    """Generate a batch of random frames and emitter positions.

    Designed exclusively for KAN/CNN-KAN training.
    Emitters are NOT sorted — hungarian_loss() handles permutation
    invariance during training so no canonical ordering is needed.

    Parameters
    ----------
    batch_size : int
    M          : number of emitters per frame
    camera     : Camera
    psf        : Gaussian2DPSF
    noise      : Noise
    Im0        : emitter intensity in photons/s
    margin_px  : pixels from each edge excluded when sampling positions

    Returns
    -------
    x    : (B, Kx*Ky) float32 — normalised flattened frames
    y_nm : (B, M, 2)  float32 — true emitter positions in nm (unsorted)
    """
    xy     = sample_emitters_batch(batch_size, M, camera, margin_px=margin_px)
    device = xy.device
    dtype  = torch.float32

    Im_tensor = torch.full(
        (batch_size, M), float(Im0), dtype=dtype, device=device
    )
    emitter = EmitterData(xy=xy, Im=Im_tensor)
    system  = SystemConfig(
        emitter = emitter,
        psf     = psf,
        camera  = camera,
        noise   = noise,
    )

    with torch.no_grad():
        frames = gauss2d_frame_torch(system)        # (B, Ky, Kx)

    x = frames.view(batch_size, -1).to(dtype)       # (B, Kx*Ky)
    x = x / (x.amax(dim=1, keepdim=True) + 1e-8)

    return x, xy   # xy: (B, M, 2) — unsorted


# ======================================================================
# Permutation-invariant Hungarian loss
# ======================================================================

def hungarian_loss(
    pred: torch.Tensor,
    true: torch.Tensor,
) -> torch.Tensor:
    """Permutation-invariant MSE loss using optimal 1-to-1 matching.

    For each sample in the batch, finds the optimal bijective assignment
    between predicted and true emitter positions (minimising total squared
    distance), then computes MSE on the matched pairs.

    Enforcing 1-to-1 assignment prevents the local-minimum problem of
    RMSMD loss where multiple predictions cluster near one true location
    while others are left unmatched — critical for M >= 3.

    Solver selected automatically (fastest available, no code change needed):

    1. torch_linear_assignment (GPU-native, any M):
       Install: pip install torch-linear-assignment
       Requires Microsoft C++ Build Tools on Windows.

    2. Brute-force permutation (pure PyTorch GPU, M <= 7):
       Enumerates all M! permutations as a pre-computed tensor.
       M=5: 120 permutations — trivially fast on GPU.
       M=7: 5040 permutations — still fast on GPU.
       No installation or compilation required.

    3. scipy fallback (CPU, any M):
       Used when M > 7 and torch_linear_assignment is unavailable.
       O(M^3) in optimised C — fine for M <= ~500.

    Parameters
    ----------
    pred : (B, M, 2) tensor — predicted positions in nm
    true : (B, M, 2) tensor — true positions in nm

    Returns
    -------
    loss : scalar tensor — mean squared error over matched pairs (nm^2)
    """
    B, M, _ = pred.shape
    device  = pred.device

    # Pairwise squared distances — (B, M, M)
    # C[b, i, j] = ||pred[b,i] - true[b,j]||^2
    C = (pred.unsqueeze(2) - true.unsqueeze(1)).pow(2).sum(dim=-1)

    if _SOLVER == 'gpu':
        # ----------------------------------------------------------
        # Branch 1: GPU-native torch_linear_assignment
        # col_ind[b,i] = j  means pred[b,i] matched to true[b,j]
        # ----------------------------------------------------------
        col_ind  = _gpu_lap(C)                                    # (B, M)
        row_idx  = torch.arange(B, device=device).unsqueeze(1).expand(B, M)
        src_idx  = torch.arange(M, device=device).unsqueeze(0).expand(B, M)
        matched_pred = pred[row_idx, src_idx]                     # (B, M, 2)
        matched_true = true[row_idx, col_ind]                     # (B, M, 2)

    elif M <= 7:
        # ----------------------------------------------------------
        # Branch 2: brute-force over all M! permutations on GPU.
        # perms : (P, M)  where P = M!
        # For each sample b and permutation p, total cost =
        #   sum_i C[b, i, perms[p, i]]
        # ----------------------------------------------------------
        perms = _get_perms(M, device)                             # (P, M)
        P     = perms.shape[0]

        # C_perm[b, p, i] = C[b, i, perms[p, i]]
        # Expand C: (B, 1, M, M) and perms: (1, P, M) for gather
        C_exp   = C.unsqueeze(1).expand(B, P, M, M)              # (B, P, M, M)
        p_exp   = perms.unsqueeze(0).unsqueeze(2).expand(B, P, 1, M)  # (B, P, 1, M)
        # We want C[b, p, i, perms[p,i]]: use diagonal trick
        # cost[b, p] = sum_i C[b, i, perms[p, i]]
        perm_exp = perms.unsqueeze(0).expand(B, P, M)            # (B, P, M)
        i_idx    = torch.arange(M, device=device).view(1, 1, M).expand(B, P, M)
        b_idx    = torch.arange(B, device=device).view(B, 1, 1).expand(B, P, M)
        cost_all = C[b_idx, i_idx, perm_exp].sum(dim=2)          # (B, P)

        # Best permutation per sample
        best_p   = cost_all.argmin(dim=1)                        # (B,)

        # Gather matched true indices using best permutation
        best_perm = perms[best_p]                                 # (B, M)
        b_idx2    = torch.arange(B, device=device).unsqueeze(1).expand(B, M)
        i_idx2    = torch.arange(M, device=device).unsqueeze(0).expand(B, M)
        matched_pred = pred[b_idx2, i_idx2]                      # (B, M, 2)
        matched_true = true[b_idx2, best_perm]                   # (B, M, 2)

    else:
        # ----------------------------------------------------------
        # Branch 3: scipy CPU fallback for M > 7
        # ----------------------------------------------------------
        C_np     = C.detach().cpu().numpy()
        matched_pred = []
        matched_true = []
        for b in range(B):
            row_ind, col_ind = _scipy_lap(C_np[b])
            matched_pred.append(pred[b][row_ind])
            matched_true.append(true[b][col_ind])
        matched_pred = torch.stack(matched_pred, dim=0)          # (B, M, 2)
        matched_true = torch.stack(matched_true, dim=0)          # (B, M, 2)

    return F.mse_loss(matched_pred, matched_true)
