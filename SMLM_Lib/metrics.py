"""
Metrics for evaluating localization performance.
"""

import torch

"""
Root Mean Square Minimum Distance (RMSMD) — a universal quality metric
for measuring the mutual fitness between two sets of points.

Reference
---------
Sun, Y. "Root Mean Square Minimum Distance as a Quality Metric for
Stochastic Optical Localization Nanoscopy Images."
Scientific Reports 8, 17211 (2018).  https://doi.org/10.1038/s41598-018-35053-8

MATLAB original: RMSMD.m  (Yi Sun, 04/08/2017, rev. 08/19/2019)
Python translation: SMLM_Lib (2025)
"""


def rmsmd(
    S: torch.Tensor,
    X: torch.Tensor,
) -> tuple:
    """Root Mean Square Minimum Distance (RMSMD) between point sets S and X.

    Computes

        D²(X, S) = [Σ_{s∈S} min_{x∈X} ‖x − s‖²
                  + Σ_{x∈X} min_{s∈S} ‖s − x‖²] / (|S| + |X|)

        D(X, S)  = √D²(X, S)

    Both summation directions are required; using only one direction can give
    misleading results when the two sets have very different cardinalities
    (see Sun 2018, Eq. 1 and surrounding discussion).

    Parameters
    ----------
    S : (M, d) tensor
        M points of dimension d, e.g., ground-truth emitter locations in nm.
    X : (N, d) tensor
        N points of dimension d, e.g., estimated emitter locations in nm.
        S and X must share the same device and dtype.

    Returns
    -------
    D  : scalar tensor — RMSMD, in the same units as S and X.
    D2 : scalar tensor — MSMD (= D²).

    Raises
    ------
    ValueError
        If either set is empty or if S and X have different point dimensions.

    Notes
    -----
    Convention change from the MATLAB original
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The MATLAB function uses shape ``(d, N)`` — each *column* is a point.
    This function uses shape ``(N, d)`` — each *row* is a point, which is
    the standard NumPy / PyTorch convention.  Convert MATLAB-style arrays
    with ``S_py = torch.tensor(S_matlab.T)``.

    Performance
    ~~~~~~~~~~~
    The full pairwise squared-distance matrix of shape ``(N, M)`` is built
    by ``torch.cdist``.  Peak memory scales as O(N · M).  For very large
    point sets (N, M ≳ 50 000 on a typical GPU) consider calling this
    function with smaller batches or moving tensors to CPU.

    Examples
    --------
    >>> import torch
    >>> S = torch.tensor([[0., 0.], [100., 0.]])      # 2 ground-truth points
    >>> X = torch.tensor([[10., 0.], [90., 0.]])      # 2 estimated points
    >>> D, D2 = rmsmd(S, X)
    >>> print(f"RMSMD = {D.item():.2f} nm")
    RMSMD = 10.00 nm
    """

    # ------------------------------------------------------------------
    # Input validation
    # ------------------------------------------------------------------
    if S.ndim != 2 or X.ndim != 2:
        raise ValueError(
            "S and X must be 2-D tensors of shape (N, d); "
            f"got S.shape={tuple(S.shape)}, X.shape={tuple(X.shape)}."
        )

    d_S = S.shape[1]
    d_X = X.shape[1]
    if d_S != d_X:
        raise ValueError(
            f"S and X must have the same point dimension d; "
            f"got S: d={d_S}, X: d={d_X}."
        )

    M = S.shape[0]
    N = X.shape[0]
    if M == 0 or N == 0:
        raise ValueError("S and X must both be non-empty.")

    # ------------------------------------------------------------------
    # Pairwise squared Euclidean distances
    # dist2[i, j] = ||X[i] - S[j]||^2,  shape (N, M)
    # ------------------------------------------------------------------
    dist2 = torch.cdist(X.float(), S.float(), p=2).pow(2)   # (N, M)

    # For each x in X: minimum squared distance to any point in S
    min_X2S = dist2.min(dim=1).values    # (N,)

    # For each s in S: minimum squared distance to any point in X
    min_S2X = dist2.min(dim=0).values    # (M,)

    # ------------------------------------------------------------------
    # MSMD and RMSMD  (Sun 2018, Eq. 1)
    # ------------------------------------------------------------------
    D2 = (min_X2S.sum() + min_S2X.sum()) / (N + M)
    D  = torch.sqrt(D2)

    return D, D2


def rmsmd_per_frame(
    pred: torch.Tensor,
    true: torch.Tensor,
) -> tuple:
    """RMSMD computed independently for each frame in a batch.

    Equivalent in semantics to calling ``rmsmd(pred[b], true[b])`` for
    every frame index ``b``, but fully vectorised over the batch dimension
    so no Python loop is needed.

    Parameters
    ----------
    pred : (B, N, d) tensor
        B frames, each with N predicted d-dimensional locations (nm).
    true : (B, M, d) tensor
        B frames, each with M ground-truth d-dimensional locations (nm).
        B and d must match those of ``pred``; N and M may differ.

    Returns
    -------
    D  : (B,) tensor — per-frame RMSMD in the same units as the inputs.
    D2 : (B,) tensor — per-frame MSMD (= D²).

    Raises
    ------
    ValueError
        If tensors are not 3-D, batch or point dimensions are inconsistent,
        or either set is empty.

    Notes
    -----
    For a single frame call ``rmsmd(pred[0], true[0])`` directly
    (see :func:`rmsmd` in this module).

    The denominator in the MSMD formula is ``N + M`` (Sun 2018, Eq. 1),
    which equals ``2*M`` only when the two sets have the same cardinality.
    """
    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------
    if pred.ndim != 3 or true.ndim != 3:
        raise ValueError(
            "pred and true must be 3-D tensors of shape (B, N, d); "
            f"got pred.shape={tuple(pred.shape)}, true.shape={tuple(true.shape)}."
        )

    B_pred, N, d_pred = pred.shape
    B_true, M, d_true = true.shape

    if B_pred != B_true:
        raise ValueError(
            f"Batch dimension B must match; got pred B={B_pred}, true B={B_true}."
        )
    if d_pred != d_true:
        raise ValueError(
            f"Point dimension d must match; got pred d={d_pred}, true d={d_true}."
        )
    if N == 0 or M == 0:
        raise ValueError("pred and true must each contain at least one point per frame.")

    # ------------------------------------------------------------------
    # Pairwise squared Euclidean distances
    # d2[b, i, j] = ||pred[b, i] - true[b, j]||^2,  shape (B, N, M)
    # ------------------------------------------------------------------
    d2 = (pred.unsqueeze(2) - true.unsqueeze(1)).pow(2).sum(dim=-1)

    # For each predicted point: minimum squared distance to any true point
    min_pred2true = d2.min(dim=2).values    # (B, N)

    # For each true point: minimum squared distance to any predicted point
    min_true2pred = d2.min(dim=1).values    # (B, M)

    # ------------------------------------------------------------------
    # Per-frame MSMD and RMSMD  (Sun 2018, Eq. 1)
    # Denominator is N + M, not 2*M — the two sets may have different sizes.
    # ------------------------------------------------------------------
    D2 = (min_pred2true.sum(dim=1) + min_true2pred.sum(dim=1)) / (N + M)
    D  = torch.sqrt(D2)

    return D, D2

# ===========================================================================
# Partition algorithm and partition-based metrics
# (Sun, Y. JOSA A, vol. 39, no. 12, pp. 2307-2315, Dec. 2022)
# MATLAB originals: partitionX.m, RMSE_P.m, RMSMD_P.m (Yi Sun, 2021-2022)
# Python translation: SMLM_Lib (2025)
# ===========================================================================


def partition_x(
    S:   torch.Tensor,
    X:   torch.Tensor,
    Kai: torch.Tensor,
) -> tuple:
    """Partition estimated locations X according to true locations S.

    Implements Algorithm 1 from Sun (2022) JOSA A.  Uses three conditions:
    (i) true locations S are known; (ii) number of activations Kai per
    emitter is known (determines prior beta_i); (iii) an estimated location
    is most likely associated with its nearest true emitter.

    For JML algorithms where K < 2*M (one estimate per emitter per frame,
    our case), equal priors beta_i = 1/M are used and Ki = 1 for all i.

    Parameters
    ----------
    S   : (M, d) tensor — M true emitter locations
    X   : (K, d) tensor — K estimated locations
    Kai : (M,)   tensor — number of activations per emitter (int)
                          For JML (K < 2M) this is ignored; equal prior used.

    Returns
    -------
    Xp  : (Kp, d) tensor — partitioned estimates, rows 0..Ki[0]-1 belong to
                           emitter 0, rows Ki[0]..Ki[0]+Ki[1]-1 to emitter 1,
                           etc.  Kp >= K (may exceed K if an emitter gets the
                           nearest estimate assigned twice in the fallback step)
    Ki  : (M,) int tensor — number of estimates assigned to each emitter
    Ip  : (Kp, 2) int tensor — Ip[n, 0] = original index in X,
                               Ip[n, 1] = emitter index assigned to

    References
    ----------
    Sun, Y. "Partition of estimated locations: an approach to accurate
    quality metrics for stochastic optical localization nanoscopy."
    JOSA A 39, 2307-2315 (2022).  https://doi.org/10.1364/JOSAA.474218
    """
    M, d = S.shape
    K    = X.shape[0]
    device = S.device
    dtype  = S.dtype

    if M == 0 or K == 0:
        raise ValueError("S and X must both be non-empty.")

    # Special case: M == 1 — all estimates belong to the single emitter
    if M == 1:
        Ki = Kai.clone()
        Xp = X.clone()
        Ip = torch.stack([
            torch.arange(K, device=device),
            torch.zeros(K, dtype=torch.long, device=device),
        ], dim=1)    # (K, 2)
        return Xp, Ki, Ip

    # ------------------------------------------------------------------
    # Pairwise squared distances  D[m, k] = ||X[k] - S[m]||^2
    # Shape: (M, K)
    # ------------------------------------------------------------------
    # X: (K, d), S: (M, d)  →  dist2: (M, K)
    D  = torch.cdist(S.float(), X.float(), p=2).pow(2)   # (M, K)
    Db = D.clone()   # backup for fallback step

    # ------------------------------------------------------------------
    # Prior and quota Ki
    # ------------------------------------------------------------------
    if K < 2 * M:                              # JML algorithm
        beta = torch.ones(M, device=device, dtype=torch.float32) / M
    else:                                      # FFL algorithm
        Ka   = Kai.float().sum()
        beta = Kai.float() / Ka

    Ki_quota = torch.round(K * beta).long()    # (M,) — quota per emitter

    # Working arrays (Python lists, converted to tensors at end)
    # Xt[m] holds the assigned estimates for emitter m
    Xt = [[] for _ in range(M)]   # list of lists of row indices into X
    pt = torch.zeros(M, dtype=torch.long, device=device)   # counter per emitter
    It_x   = []   # original X index
    It_emitter = []   # emitter index

    TKi = Ki_quota.sum().item()

    D = D.clone().float()   # working copy (will be modified with Inf)
    inf = float('inf')

    def _find_min_pair(D_mat):
        """Find (m*, k*) = argmin_{m,k} D_mat[m,k]."""
        flat_idx = torch.argmin(D_mat)
        m = (flat_idx // K).item()
        c = (flat_idx %  K).item()
        return m, c

    n = 0   # assignment counter

    if K <= TKi:
        # ------------------------------------------------------------------
        # Case K <= TKi: assign all K estimates (Algorithm 1, step ii)
        # ------------------------------------------------------------------
        for _ in range(K):
            m, c = _find_min_pair(D)
            pt[m] += 1
            Xt[m].append(c)
            D[:, c] = inf          # remove estimate c
            if pt[m] == Ki_quota[m]:
                D[m, :] = inf      # emitter m's quota filled
            It_x.append(c)
            It_emitter.append(m)
            n += 1
    else:
        # ------------------------------------------------------------------
        # Case K > TKi: two-pass assignment (Algorithm 1, steps iii-iv)
        # Pass 1: assign TKi estimates respecting quota
        # ------------------------------------------------------------------
        Dt = D.clone()
        for _ in range(int(TKi)):
            m, c = _find_min_pair(D)
            pt[m] += 1
            Xt[m].append(c)
            D[:, c]  = inf
            if pt[m] == Ki_quota[m]:
                D[m, :] = inf
            Dt[:, c] = inf         # remove from Dt too (but keep emitter rows)
            It_x.append(c)
            It_emitter.append(m)
            n += 1

        # Pass 2: assign remaining K - TKi estimates (one per emitter, no quota)
        D = Dt
        for _ in range(K - int(TKi)):
            m, c = _find_min_pair(D)
            pt[m] += 1
            Xt[m].append(c)
            D[:, c] = inf
            D[m, :] = inf          # remove emitter too (K-TKi <= M)
            It_x.append(c)
            It_emitter.append(m)
            n += 1

    # ------------------------------------------------------------------
    # Fallback: if emitter m got no estimate, assign its nearest (step v)
    # ------------------------------------------------------------------
    for m in range(M):
        if pt[m] == 0:
            c = int(torch.argmin(Db[m, :]).item())
            pt[m] += 1
            Xt[m].append(c)
            It_x.append(c)
            It_emitter.append(m)
            n += 1

    Ki_out = pt   # (M,) — actual counts after partition

    # ------------------------------------------------------------------
    # Build Xp — partitioned estimates stacked by emitter
    # ------------------------------------------------------------------
    Xp_rows = []
    for m in range(M):
        for c in Xt[m]:
            Xp_rows.append(X[c])
    Xp = torch.stack(Xp_rows, dim=0)   # (Kp, d)

    Ip = torch.stack([
        torch.tensor(It_x[:n],       dtype=torch.long, device=device),
        torch.tensor(It_emitter[:n], dtype=torch.long, device=device),
    ], dim=1)   # (Kp, 2)

    return Xp, Ki_out, Ip


def rmse_p(
    S:  torch.Tensor,
    Xp: torch.Tensor,
    Ki: torch.Tensor,
) -> torch.Tensor:
    """Sample RMSE with known partition (RMSE-P).

    Computes the sample root mean square error between estimated locations
    and their assigned true emitters, using the partition produced by
    :func:`partition_x`.

        h²(X, S) = (1/K) Σ_m Σ_{x ∈ X_m} ||x − s_m||²

        h(X, S)  = √h²(X, S)

    Parameters
    ----------
    S  : (M, d) tensor — true emitter locations
    Xp : (K, d) tensor — partitioned estimates (output of partition_x)
         rows 0..Ki[0]-1 belong to emitter 0, etc.
    Ki : (M,) int tensor — number of estimates per emitter (from partition_x)

    Returns
    -------
    h  : scalar tensor — RMSE-P in the same units as S and Xp

    References
    ----------
    Sun, Y. JOSA A 39, 2307-2315 (2022), Eq. (3) with estimated partition.
    """
    M = S.shape[0]
    K = Xp.shape[0]

    h2 = torch.tensor(0.0, device=S.device, dtype=torch.float32)
    p  = 0
    for m in range(M):
        k_m = int(Ki[m].item())
        if k_m > 0:
            diff = Xp[p : p + k_m].float() - S[m].float().unsqueeze(0)  # (k_m, d)
            h2   = h2 + diff.pow(2).sum()
        p += k_m

    return torch.sqrt(h2 / K)


def rmsmd_p(
    S:  torch.Tensor,
    Xp: torch.Tensor,
    Ki: torch.Tensor,
) -> tuple:
    """RMSMD with known partition (RMSMD-P).

    Computes the modified RMSMD using the partition produced by
    :func:`partition_x` (Sun 2022, Eq. 11):

        Dp²(X, S) = (1/(K+M)) Σ_m [ min_{x ∈ X_m} ||x − s_m||²
                                    + Σ_{x ∈ X_m} ||x − s_m||² ]

    For JML algorithms (K = M, one estimate per emitter), this equals
    RMSE-P exactly (Sun 2022, Eq. 14).

    Parameters
    ----------
    S  : (M, d) tensor — true emitter locations
    Xp : (K, d) tensor — partitioned estimates (output of partition_x)
    Ki : (M,) int tensor — number of estimates per emitter

    Returns
    -------
    D  : scalar tensor — RMSMD-P in the same units as S and Xp
    D2 : scalar tensor — MSMD-P (= D²)

    References
    ----------
    Sun, Y. JOSA A 39, 2307-2315 (2022), Eq. (11).
    """
    M = S.shape[0]
    K = Xp.shape[0]

    D_X2S = torch.tensor(0.0, device=S.device, dtype=torch.float32)
    D_S2X = torch.tensor(0.0, device=S.device, dtype=torch.float32)
    p = 0
    for m in range(M):
        k_m = int(Ki[m].item())
        if k_m > 0:
            diff = Xp[p : p + k_m].float() - S[m].float().unsqueeze(0)  # (k_m, d)
            d2_m = diff.pow(2).sum(dim=1)   # (k_m,) — per-estimate sq dist to s_m
            D_X2S = D_X2S + d2_m.sum()     # Σ_{x ∈ X_m} ||x - s_m||^2
            D_S2X = D_S2X + d2_m.min()     # min_{x ∈ X_m} ||x - s_m||^2
        p += k_m

    D2 = (D_S2X + D_X2S) / (K + M)
    D  = torch.sqrt(D2)
    return D, D2
