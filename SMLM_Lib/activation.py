"""
activation.py
Emitter photoactivation processes for SMLM data movies.

Implements the Markov-chain models of emitter activation from Sun (2024,
Optics Express, "Markov chain models of emitter activations in single
molecule localization microscopy") for both continuous and cycled
illuminations. These are the activation sub-model of the data-movie model;
combined with the data-frame model (e.g. ``gauss2d_frame_torch``) they
produce a full data movie.

State convention (per frame, per emitter), following the paper:
    0 - deactivated   (photoactivatable, no photon emission)
    1 - activated     (photoactivatable, emitting photons)
    2 - photobleached (absorbing state)

Both functions return the state matrix ``ca`` of shape ``(N, M)`` -- frame
index by emitter index -- following the (N, M) tensor convention used
throughout SMLM_Lib (see ``EmitterData``). To drive the frame model, build
the per-frame intensity matrix directly::

    Im = torch.where(ca == 1, Im0, 0.0)   # (N, M); intensity on where state==1

with the M emitters held at fixed structure positions broadcast across the
N frames.

These routines are torch ports of:
    emActMarkovContinue.m  (04/29/2024)
    emActMarkovCycle.m     (04/26/2024)
adapted to the (N, M) array convention of SMLM_Lib.

Notes
-----
- All times are in seconds.
- The simulation is sequential over frames N (each frame depends on the
  previous) and fully vectorised over the M emitters.
- An emitter in the absorbing state 2 needs no special-casing: gathering
  column 2 of the transition matrix, [0, 0, 1]^T, makes the inverse-CDF
  step return state 2 with probability 1.
- A single uniform initial draw assigns each emitter to state 0 or 1
  according to the stationary distribution, then one transition is applied
  before frame 1 is recorded (the burn-in row ``c[0]`` is dropped).
"""

import math
from typing import Optional, Tuple

import torch


# ============================================================
# Internal helper: one vectorised Markov transition step
# ============================================================

def _step(cur: torch.Tensor, R: torch.Tensor, gen) -> torch.Tensor:
    """Advance M emitters one frame under transition matrix ``R``.

    Parameters
    ----------
    cur : (M,) long tensor
        Current states in {0, 1, 2}.
    R : (3, 3) tensor
        Column-stochastic transition matrix: ``R[j, i]`` is the probability
        of moving from state i to state j.
    gen : torch.Generator or None
        RNG for reproducibility.

    Returns
    -------
    nxt : (M,) long tensor
        Next states in {0, 1, 2}, sampled by inverse-CDF on the column of
        ``R`` selected by each emitter's current state.
    """
    cols = R[:, cur]                                   # (3, M): per-emitter column
    thr0 = cols[0]                                     # P(next == 0)
    thr1 = cols[0] + cols[1]                           # P(next in {0, 1})
    P = torch.rand(cur.shape[0], generator=gen,
                   dtype=R.dtype, device=R.device)     # (M,) ~ U[0, 1)
    # next = 0 if P < thr0; 1 if thr0 <= P < thr1; 2 if P >= thr1.
    # For state 2 the column is [0, 0, 1] -> thr0 = thr1 = 0 -> next = 2.
    return (P >= thr0).long() + (P >= thr1).long()


def _resolve_device(device) -> torch.device:
    if device is None:
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    return torch.device(device)


def _cycle_stationary_h1(A: torch.Tensor, D: torch.Tensor,
                         KA: int, KD: int) -> torch.Tensor:
    """Per-frame stationary activated-state probability for a two-phase cycle.

    Computes h^(k)_1 for k = 1, ..., K (Sun 2024, Eqs. 13-14) as t -> inf,
    i.e. ignoring photobleaching. For the k-th frame of a cycle the one-cycle
    transition matrix Q^(k) is assembled from the 2x2 phase matrices A and D,
    and the stationary activated probability is

        h^(k)_1 = (1 - Q^(k)_00) / (2 - Q^(k)_00 - Q^(k)_11).

    Parameters
    ----------
    A, D : (2, 2) column-stochastic phase matrices (A-phase, D-phase).
    KA, KD : frames per phase (K = KA + KD).

    Returns
    -------
    h1 : (K,) tensor — h^(k)_1 indexed by k-1 (k = 1..K).
    """
    K = KA + KD
    mp = torch.linalg.matrix_power
    h1 = torch.zeros(K, dtype=A.dtype, device=A.device)
    for i in range(1, K + 1):
        if i <= KA:
            Q = mp(A, i) @ mp(D, KD) @ mp(A, KA - i)
        else:
            Q = mp(D, i - KA) @ mp(A, KA) @ mp(D, K - i)
        h1[i - 1] = (1.0 - Q[0, 0]) / (2.0 - Q[0, 0] - Q[1, 1])
    return h1


# ============================================================
# Continuous illumination
# ============================================================

def emActMarkovContinue(
    t: float,
    t1: float,
    t0: float,
    Dt: float,
    N: int,
    M: int,
    device=None,
    seed: Optional[int] = None,
) -> Tuple[torch.Tensor, float, float]:
    """Emitter activation states for a continuous-illumination data movie.

    Single-phase three-state Markov chain (Sun 2024, Sec. 3). The same
    transition matrix is applied in every frame.

    Parameters
    ----------
    t  : mean of the photoactivatable period (s).
    t1 : mean of the activation period (s), t1 << t.
    t0 : mean of the deactivation period (s), t0 << t, t0 > t1.
    Dt : frame time (s), Dt << t0, t.
    N  : number of frames in the data movie.
    M  : number of emitters.
    device : torch device or None (auto: CUDA if available, else CPU).
    seed   : optional int for reproducible RNG.

    Returns
    -------
    ca : (N, M) long tensor
        State of emitter m in frame n; values in {0, 1, 2}.
    p  : float
        Per-frame probability of remaining photoactivatable, exp(-Dt/t).
    h1 : float
        Stationary probability of the activated state as t -> inf.
    """
    device = _resolve_device(device)
    gen = (torch.Generator(device=device).manual_seed(int(seed))
           if seed is not None else None)

    p0 = math.exp(-Dt / t0)        # retain state 0 over one frame
    p1 = math.exp(-Dt / t1)        # retain state 1 over one frame
    p = math.exp(-Dt / t)          # remain photoactivatable over one frame
    h1 = (1.0 - p0) / (2.0 - p0 - p1)

    # Transition probabilities (with bleaching factor p).
    r00 = p0 * p
    r10 = (1.0 - p0) * p
    r20 = 1.0 - p
    r01 = (1.0 - p1) * p
    r11 = p1 * p
    r21 = 1.0 - p
    R = torch.tensor(
        [[r00, r01, 0.0],
         [r10, r11, 0.0],
         [r20, r21, 1.0]],
        dtype=torch.float64, device=device,
    )

    # Row c[0] is the stationary-distribution initial draw (burn-in).
    c = torch.zeros(N + 1, M, dtype=torch.long, device=device)
    u0 = torch.rand(M, generator=gen, dtype=torch.float64, device=device)
    c[0] = (u0 <= h1).long()                          # state 1 w.p. h1, else 0

    for n in range(N):
        c[n + 1] = _step(c[n], R, gen)

    ca = c[1:N + 1].contiguous()                       # (N, M); drop burn-in
    return ca, float(p), float(h1)


# ============================================================
# Cycled illumination (two-phase)
# ============================================================

def emActMarkovCycle(
    t: float,
    tA1: float,
    tA0: float,
    tD1: float,
    tD0: float,
    Dt: float,
    C: int,
    KA: int,
    KD: int,
    M: int,
    device=None,
    seed: Optional[int] = None,
) -> Tuple[torch.Tensor, float, torch.Tensor]:
    """Emitter activation states for a cycled-illumination data movie.

    Two-phase three-state Markov chain (Sun 2024, Sec. 2). Each cycle is an
    A-phase of KA activation frames followed by a D-phase of KD deactivation
    frames; the movie is C cycles, so N = C * (KA + KD) frames.

    Parameters
    ----------
    t   : mean of the photoactivatable period (s).
    tA1 : mean of the activation period in the A-phase (s), tA1 << t.
    tA0 : mean of the deactivation period in the A-phase (s), tA0 << t,
          tA0 > tA1.
    tD1 : mean of the activation period in the D-phase (s), tD1 << t.
    tD0 : mean of the deactivation period in the D-phase (s), tD0 << t,
          tD0 > tD1. Typically tA0 << tD0 and tA1 >> tD1.
    Dt  : frame time (s), Dt << tA0, tD0, t.
    C   : number of cycles.
    KA  : number of A-phase frames per cycle.
    KD  : number of D-phase frames per cycle (K = KA + KD, N = C * K).
    M   : number of emitters.
    device : torch device or None (auto: CUDA if available, else CPU).
    seed   : optional int for reproducible RNG.

    Returns
    -------
    ca : (N, M) long tensor
        State of emitter m in frame n; values in {0, 1, 2}.
    p  : float
        Per-frame probability of remaining photoactivatable, exp(-Dt/t).
    h1 : (K,) tensor
        Stationary probability of the activated state in each of the K
        frames of a cycle as t -> inf.
    """
    device = _resolve_device(device)
    gen = (torch.Generator(device=device).manual_seed(int(seed))
           if seed is not None else None)

    K = KA + KD
    N = C * K

    pA0 = math.exp(-Dt / tA0)
    pA1 = math.exp(-Dt / tA1)
    pD0 = math.exp(-Dt / tD0)
    pD1 = math.exp(-Dt / tD1)
    p = math.exp(-Dt / t)

    # Phase matrices as t -> inf (2x2, column-stochastic) used for the
    # per-cycle stationary distribution of the activated state.
    A = torch.tensor([[pA0, 1.0 - pA1],
                      [1.0 - pA0, pA1]], dtype=torch.float64, device=device)
    D = torch.tensor([[pD0, 1.0 - pD1],
                      [1.0 - pD0, pD1]], dtype=torch.float64, device=device)

    h1 = _cycle_stationary_h1(A, D, KA, KD)    # (K,) per-frame stationary, Eq (14)

    # Per-phase three-state transition matrices (with bleaching factor p).
    a00, a10, a20 = pA0 * p, (1.0 - pA0) * p, 1.0 - p
    a01, a11, a21 = (1.0 - pA1) * p, pA1 * p, 1.0 - p
    d00, d10, d20 = pD0 * p, (1.0 - pD0) * p, 1.0 - p
    d01, d11, d21 = (1.0 - pD1) * p, pD1 * p, 1.0 - p
    Ra = torch.tensor([[a00, a01, 0.0],
                       [a10, a11, 0.0],
                       [a20, a21, 1.0]], dtype=torch.float64, device=device)
    Rd = torch.tensor([[d00, d01, 0.0],
                       [d10, d11, 0.0],
                       [d20, d21, 1.0]], dtype=torch.float64, device=device)

    c = torch.zeros(N + 1, M, dtype=torch.long, device=device)
    u0 = torch.rand(M, generator=gen, dtype=torch.float64, device=device)
    c[0] = (u0 <= h1[K - 1]).long()        # init from end-of-cycle stationary

    n = 0
    for _ in range(C):
        for _i in range(KA):               # A-phase
            c[n + 1] = _step(c[n], Ra, gen)
            n += 1
        for _i in range(KD):               # D-phase
            c[n + 1] = _step(c[n], Rd, gen)
            n += 1

    ca = c[1:N + 1].contiguous()           # (N, M); drop burn-in
    return ca, float(p), h1


# ============================================================
# Analytical predictions (no simulation)
# ============================================================

def activation_statistics_continue(
    t: float,
    t1: float,
    t0: float,
    Dt: float,
    N: int,
    M: int,
) -> dict:
    """Analytical activation statistics for continuous illumination.

    Closed-form predictions from Sun (2024), Sec. 3.3, to compare against
    a simulation produced by :func:`emActMarkovContinue` with the same
    parameters. All quantities account for photobleaching (p < 1).

    Parameters match :func:`emActMarkovContinue` (t, t1, t0, Dt, N, M).

    Returns
    -------
    dict with keys:
        'p', 'p0', 'p1' : per-frame dwelling probabilities.
        'h1'  : stationary activated-state probability (t -> inf), Eq. (25).
        'Np'  : avg photoactivatable frames per emitter, Eq. (15).
        'Nae' : avg activated frames per emitter, Eq. (28).
        'Na'  : avg total activated frames over all M emitters, Eq. (29).
        'Ma'  : (N,) tensor, avg activated emitters in each frame, Eq. (27).
    """
    p0 = math.exp(-Dt / t0)
    p1 = math.exp(-Dt / t1)
    p = math.exp(-Dt / t)
    h1 = (1.0 - p0) / (2.0 - p0 - p1)

    n = torch.arange(1, N + 1, dtype=torch.float64)
    Ma = M * (p ** n) * h1                                  # Eq (27)

    if abs(1.0 - p) < 1e-12:                                # no bleaching limit
        Np = float(N)
    else:
        Np = p * (1.0 - p ** N) / (1.0 - p)                 # Eq (15)
    Nae = Np * h1                                           # Eq (28)
    Na = M * Nae                                            # Eq (29)

    return {'p': p, 'p0': p0, 'p1': p1, 'h1': h1,
            'Np': Np, 'Nae': Nae, 'Na': Na, 'Ma': Ma}


def activation_statistics_cycle(
    t: float,
    tA1: float,
    tA0: float,
    tD1: float,
    tD0: float,
    Dt: float,
    C: int,
    KA: int,
    KD: int,
    M: int,
) -> dict:
    """Analytical activation statistics for cycled illumination.

    Closed-form predictions from Sun (2024), Secs. 2.2-2.4, to compare
    against a simulation produced by :func:`emActMarkovCycle` with the same
    parameters. All quantities account for photobleaching (p < 1).

    Parameters match :func:`emActMarkovCycle`
    (t, tA1, tA0, tD1, tD0, Dt, C, KA, KD, M).

    Returns
    -------
    dict with keys:
        'p'    : per-frame photoactivatable dwelling probability.
        'h1'   : (K,) tensor, per-cycle-frame stationary activated
                 probability h^(k)_1 (t -> inf), Eq. (14).
        'Np'   : avg photoactivatable frames per emitter, Eq. (15).
        'Nae1' : avg activated frames per emitter in the first cycle, Eq. (19).
        'Nae'  : avg activated frames per emitter over the movie, Eq. (18).
        'Na'   : avg total activated frames over all M emitters, Eq. (20).
        'Ma'   : (N,) tensor, avg activated emitters in each frame, Eq. (17).
    """
    K = KA + KD
    N = C * K

    pA0 = math.exp(-Dt / tA0)
    pA1 = math.exp(-Dt / tA1)
    pD0 = math.exp(-Dt / tD0)
    pD1 = math.exp(-Dt / tD1)
    p = math.exp(-Dt / t)

    A = torch.tensor([[pA0, 1.0 - pA1],
                      [1.0 - pA0, pA1]], dtype=torch.float64)
    D = torch.tensor([[pD0, 1.0 - pD1],
                      [1.0 - pD0, pD1]], dtype=torch.float64)
    h1 = _cycle_stationary_h1(A, D, KA, KD)                 # (K,), Eq (14)

    # Ma(n) for n = (l-1)K + k: position k within the cycle picks h^(k)_1.
    n = torch.arange(1, N + 1, dtype=torch.float64)
    k_idx = (n.long() - 1) % K                              # 0-based, (N,)
    Ma = M * (p ** n) * h1[k_idx]                           # Eq (17)

    k = torch.arange(1, K + 1, dtype=torch.float64)
    Nae1 = torch.sum((p ** k) * h1).item()                  # Eq (19)

    if abs(1.0 - p) < 1e-12:                                # no bleaching limit
        Np = float(N)
        Nae = float(C * Nae1)                               # Nae ~= C * Nae(1)
    else:
        Np = p * (1.0 - p ** N) / (1.0 - p)                 # Eq (15)
        Nae = (1.0 - p ** N) / (1.0 - p ** K) * Nae1        # Eq (18)
    Na = M * Nae                                            # Eq (20)

    return {'p': p, 'h1': h1, 'Np': Np, 'Nae1': Nae1,
            'Nae': Nae, 'Na': Na, 'Ma': Ma}
