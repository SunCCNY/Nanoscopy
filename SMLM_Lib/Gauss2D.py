"""
Gauss2D.py
Functions related to the 2D Gaussian PSF.

Implements the data-frame model of Sun (2013, JBO) and Sun & Guan (2021,
JOSA A) with the per-pixel-noise extension and the per-frame variable
emitter-locations extension adopted by SMLM_Lib.
"""

import math
from typing import Optional

import torch

from .config import (
    EmitterData,
    Gaussian2DPSF,
    Camera,
    Noise,
    SystemConfig,
)


# ============================================================
# Helpers
# ============================================================

def Qfunc_torch(x: torch.Tensor) -> torch.Tensor:
    """Gaussian Q-function

        Q(x) = 0.5 * (1 - erf(x / sqrt(2))) = 1 - Phi(x).

    Used for pixel integration of the 2D Gaussian PSF. Equivalent to
    ``0.5 * torch.special.erfc(x / sqrt(2))``.
    """
    return 0.5 * (1.0 - torch.erf(x / math.sqrt(2.0)))


# ============================================================
# Frame simulator
# ============================================================

def gauss2d_frame_torch(
    system: SystemConfig,
    chunk_size: Optional[int] = None,
) -> torch.Tensor:
    """Generate SMLM frames using the 2D Gaussian PSF.

    Implements the per-pixel mean photon count

        v(n, ky, kx) = Dt * Dx * Dy * sum_m  Im[n, m] * q_m(n, ky, kx)
                       + Dt * Dx * Dy * b(ky, kx)

    where q_m(n, ky, kx) is the pixel-averaged 2D Gaussian PSF, evaluated
    via the standard normal CDF along x and y separately. Poisson photon
    counts are sampled from v, and Gaussian readout noise (per-pixel mean
    mu and variance G) is added.

    Parameters
    ----------
    system : SystemConfig
        Must contain a Gaussian2DPSF, a 2-D EmitterData, a Camera, and a
        Noise. The emitter positions, intensities, and noise tensors may
        live on any devices; the function aligns them to ``emitter.xy``'s
        device internally.
    chunk_size : int, optional
        If given, frames are processed ``chunk_size`` at a time along the
        N axis. Useful for keeping the (N, M, Ky*Kx) intermediate tensor
        within memory for large movies or large frames. ``None`` (default)
        processes all N frames at once.

    Returns
    -------
    U : (N, Ky, Kx) tensor
        Simulated frames in photoelectron units (float). Cast to uint16
        downstream if you want camera-style output.
    """
    emitter = system.emitter
    psf = system.psf
    camera = system.camera
    noise = system.noise

    # ------------------------------------------------------------
    # Validate PSF type
    # ------------------------------------------------------------
    if psf.psf_type != "gaussian2d":
        raise ValueError(
            f"gauss2d_frame_torch requires Gaussian2DPSF; got '{psf.psf_type}'."
        )

    # ------------------------------------------------------------
    # Unpack
    # ------------------------------------------------------------
    xy = emitter.positions          # (N, M, 2)
    Im = emitter.Im                 # (N, M)
    N, M, _ = xy.shape

    device = xy.device
    dtype = xy.dtype

    Kx, Ky = camera.Kx, camera.Ky
    Dx, Dy = camera.Dx, camera.Dy
    Dt = camera.Dt
    sigma = psf.sigma

    # ------------------------------------------------------------
    # Per-pixel noise quantities — align to the emitter's device.
    # (These are (Ky, Kx) tensors after the per-pixel-noise extension;
    # callers may have constructed them on a different device.)
    # ------------------------------------------------------------
    b_pix  = noise.poisson_mean_per_pixel(camera).to(device=device, dtype=dtype)
    g_mean = noise.gaussian_mean_per_pixel(camera).to(device=device, dtype=dtype)
    g_var  = noise.gaussian_var_per_pixel(camera).to(device=device, dtype=dtype)
    g_std  = torch.sqrt(g_var)

    # ------------------------------------------------------------
    # Pixel grid boundaries (independent of frames/emitters)
    #   shape (1, 1, 1, Kx) and (1, 1, Ky, 1) for clean broadcasting
    # ------------------------------------------------------------
    kx = torch.arange(Kx, device=device, dtype=dtype).view(1, 1, 1, Kx)
    ky = torch.arange(Ky, device=device, dtype=dtype).view(1, 1, Ky, 1)

    Dxk0 = Dx * kx
    Dxk1 = Dx * (kx + 1)
    Dyk0 = Dy * ky
    Dyk1 = Dy * (ky + 1)

    # ------------------------------------------------------------
    # Allocate output and process N in chunks (chunk_size=None -> one chunk)
    # ------------------------------------------------------------
    U = torch.empty((N, Ky, Kx), device=device, dtype=dtype)
    chunk = N if chunk_size is None else int(chunk_size)
    if chunk <= 0:
        raise ValueError("chunk_size must be a positive integer or None.")

    for n0 in range(0, N, chunk):
        n1 = min(n0 + chunk, N)
        Nc = n1 - n0

        # Emitter coordinates for this chunk: (Nc, M, 1, 1)
        x0 = xy[n0:n1, :, 0].view(Nc, M, 1, 1)
        y0 = xy[n0:n1, :, 1].view(Nc, M, 1, 1)

        # Normalised pixel-edge offsets
        Dxk0_n = (Dxk0 - x0) / sigma
        Dxk1_n = (Dxk1 - x0) / sigma
        Dyk0_n = (Dyk0 - y0) / sigma
        Dyk1_n = (Dyk1 - y0) / sigma

        # Pixel-integrated 2D Gaussian PSF: Phi(upper) - Phi(lower).
        # Using Q(-x) = Phi(x) gives the same result as torch.special.erfc.
        qx = Qfunc_torch(-Dxk1_n) - Qfunc_torch(-Dxk0_n)   # (Nc, M, 1, Kx)
        qy = Qfunc_torch(-Dyk1_n) - Qfunc_torch(-Dyk0_n)   # (Nc, M, Ky, 1)
        q  = qy * qx                                         # (Nc, M, Ky, Kx)

        # Weighted sum over emitters via bmm to avoid materialising the
        # full (Nc, M, Ky, Kx) -> sum tensor explicitly.
        q_flat   = q.reshape(Nc, M, Ky * Kx)                # (Nc, M, Ky*Kx)
        Im_chunk = Im[n0:n1].view(Nc, 1, M)                 # (Nc, 1, M)
        weighted = Dt * torch.bmm(Im_chunk, q_flat).view(Nc, Ky, Kx)

        # Mean photon count per pixel
        v = weighted + b_pix                                 # (Nc, Ky, Kx)

        # Poisson sampling, then Gaussian readout
        V = torch.poisson(v)
        U[n0:n1] = V + g_std * torch.randn_like(V) + g_mean

    return U


# ============================================================
# Fisher Information Matrix
# ============================================================

def gauss2d_fisher_torch(system: SystemConfig, ImType: str) -> torch.Tensor:
    """Compute the Fisher Information Matrix (FIM) for 2D Gaussian PSF.

    Implements Theorem 1 (Poisson noise only) and Corollary 2 (Gaussian-
    noise approximation) of Sun & Guan (2021, JOSA A), Eqs (11)–(15) and
    (28)–(30), extended to the per-pixel noise model of this library.

    The parameter ordering within F is:
        locations : (x1, y1, x2, y2, ..., xM, yM)   — indices 0 .. 2M-1
        intensity : I           for UI  — index 2M
                    I1,..,IM    for UA  — indices 2M .. 3M-1

    Parameters
    ----------
    system : SystemConfig
        Must contain a Gaussian2DPSF, a 2-D EmitterData (frame index 0 is
        used — CRLB is a per-frame quantity), a Camera, and a Noise with
        per-pixel b and G.
    ImType : str
        'KI' — Known Identical intensities (Im = I, known).
               F is (2M, 2M).
        'KA' — Known Arbitrary intensities (Im known, arbitrary).
               F is (2M, 2M).
        'UI' — Unknown, known-to-be-Identical intensity (UKI in paper).
               F is (2M+1, 2M+1).
        'UA' — Unknown Arbitrary intensities (UAA in paper).
               F is (3M, 3M).

    Returns
    -------
    F : tensor
        Fisher information matrix in SI units consistent with nm and
        photons/s (specifically nm^{-2} s^{-1} for location–location
        elements).

    Notes
    -----
    The mean mu of Gaussian noise does not appear in the FIM (Sun 2013,
    JBO, Section 2.2.3) and is not used here.

    References
    ----------
    Sun & Guan, JOSA A 38 (2021), Eqs (11)–(15), (28)–(30).
    Sun, JBO 18 (2013), Eqs (22), (35)–(38), (45)–(46).
    """
    if ImType not in ('KI', 'KA', 'UI', 'UA'):
        raise ValueError(
            f"ImType must be 'KI', 'KA', 'UI', or 'UA'; got '{ImType}'."
        )

    emitter = system.emitter
    psf     = system.psf
    camera  = system.camera
    noise   = system.noise

    if psf.psf_type != "gaussian2d":
        raise ValueError("gauss2d_fisher_torch requires Gaussian2DPSF.")

    # ----------------------------------------------------------------
    # Unpack — use frame 0 (CRLB is a per-frame quantity).
    # ----------------------------------------------------------------
    xy = emitter.positions[0]       # (M, 2) nm
    Im = emitter.Im[0]              # (M,)   photons/s
    M  = xy.shape[0]

    device = xy.device
    dtype  = xy.dtype

    sigma = psf.sigma
    Kx, Ky = camera.Kx, camera.Ky
    Dx, Dy = camera.Dx, camera.Dy
    Dt = camera.Dt

    # ----------------------------------------------------------------
    # Noise maps — align to emitter device and dtype.
    # mu is not used (it does not affect the FIM).
    # ----------------------------------------------------------------
    b = noise.b.to(device=device, dtype=dtype)   # (Ky, Kx)
    G = noise.G.to(device=device, dtype=dtype)   # (Ky, Kx)

    # ----------------------------------------------------------------
    # Intensity fractions — pre-computed once before the emitter loop.
    #   For KI / UI : all βm = 1.
    #   For KA / UA : βm = Im_m / I_mean  (2021, Eq 2).
    # ----------------------------------------------------------------
    I_mean = Im.mean()                                           # scalar tensor
    beta = (Im / I_mean) if ImType in ('KA', 'UA') \
           else torch.ones(M, device=device, dtype=dtype)        # (M,)

    # ----------------------------------------------------------------
    # Pixel grid (1-D; torch.outer handles the 2-D cross products).
    # ----------------------------------------------------------------
    kx = torch.arange(Kx, device=device, dtype=dtype)   # (Kx,)
    ky = torch.arange(Ky, device=device, dtype=dtype)   # (Ky,)

    Dxk0 = Dx * kx          # (Kx,)  left edge of pixel kx
    Dxk1 = Dx * (kx + 1)    # (Kx,)  right edge
    Dyk0 = Dy * ky           # (Ky,)
    Dyk1 = Dy * (ky + 1)    # (Ky,)

    sqrt2pi_sigma = math.sqrt(2.0 * math.pi) * sigma

    # ----------------------------------------------------------------
    # Per-emitter arrays — computed once, shared across all FIM cases.
    #
    #   q_all : (Ky, Kx, M)  pixel-averaged PSF qm(kx,ky)        [1/nm^2]
    #   DqDx  : (Ky, Kx, M)  ∂qm/∂xm                            [1/nm^3]
    #   DqDy  : (Ky, Kx, M)  ∂qm/∂ym                            [1/nm^3]
    #   Q     : (Ky, Kx)     Σm βm qm(kx,ky)                    [1/nm^2]
    # ----------------------------------------------------------------
    q_all = torch.zeros(Ky, Kx, M, device=device, dtype=dtype)
    DqDx  = torch.zeros(Ky, Kx, M, device=device, dtype=dtype)
    DqDy  = torch.zeros(Ky, Kx, M, device=device, dtype=dtype)
    Q     = torch.zeros(Ky, Kx,    device=device, dtype=dtype)

    for i in range(M):
        x0 = xy[i, 0]
        y0 = xy[i, 1]

        # Normalised pixel-edge offsets (dimensionless)
        Dxk0_n = (Dxk0 - x0) / sigma    # (Kx,)
        Dxk1_n = (Dxk1 - x0) / sigma
        Dyk0_n = (Dyk0 - y0) / sigma    # (Ky,)
        Dyk1_n = (Dyk1 - y0) / sigma

        # Pixel-averaged PSF marginals [1/nm] — Sun 2013, Eqs (45)–(46).
        # q_x(kx) = [Phi(Dxk1_n) - Phi(Dxk0_n)] / Dx, using Q(-u)=Phi(u).
        qx = (Qfunc_torch(-Dxk1_n) - Qfunc_torch(-Dxk0_n)) / Dx   # (Kx,)
        qy = (Qfunc_torch(-Dyk1_n) - Qfunc_torch(-Dyk0_n)) / Dy   # (Ky,)

        # PSF x-derivative [1/nm^2] — Sun 2013, Eq (37).
        # dq_x/dxm = (1/Dx) * [-phi(Dxk1_n) + phi(Dxk0_n)]
        DqxDx = (torch.exp(-Dxk0_n ** 2 / 2.0)
                 - torch.exp(-Dxk1_n ** 2 / 2.0)) / (Dx * sqrt2pi_sigma)
        DqyDy = (torch.exp(-Dyk0_n ** 2 / 2.0)
                 - torch.exp(-Dyk1_n ** 2 / 2.0)) / (Dy * sqrt2pi_sigma)

        # 2-D PSF and derivatives via separability — Sun 2013, Eqs (35)–(36).
        qxy  = torch.outer(qy, qx)        # (Ky, Kx)
        dqdx = torch.outer(qy, DqxDx)    # (Ky, Kx)
        dqdy = torch.outer(DqyDy, qx)    # (Ky, Kx)

        q_all[:, :, i] = qxy
        DqDx[:, :, i]  = dqdx
        DqDy[:, :, i]  = dqdy

        # Accumulate Q = Σm βm qm(k) — 2021, Eq (2).
        Q += beta[i] * qxy

    # ----------------------------------------------------------------
    # Effective denominator: Q(k) + [b(k) + G(k)] / I_mean
    #   = Q(k) + γp^{-1}(k) + γg^{-1}(k)     — 2021, Eq (28).
    # Clamped for numerical safety (both b >= 0 and G >= 0 by construction,
    # but Q can be ~0 at pixels far from all emitters).
    # ----------------------------------------------------------------
    Qu = torch.clamp(Q + (b + G) / I_mean, min=1e-10)  # (Ky, Kx)

    DtDxDy = Dt * Dx * Dy

    # ================================================================
    # CASE KI — all βm = 1, intensities known.  F = F11,  (2M × 2M).
    # Sun & Guan 2021, Eq (28) with βi = βm = 1.
    # ================================================================
    if ImType == 'KI':
        F = torch.zeros(2 * M, 2 * M, device=device, dtype=dtype)

        for i in range(M):
            for j in range(i, M):
                ii, jj = 2 * i, 2 * j

                Fxx = torch.sum(DqDx[:, :, i] * DqDx[:, :, j] / Qu)
                Fxy = torch.sum(DqDx[:, :, i] * DqDy[:, :, j] / Qu)
                Fyx = torch.sum(DqDy[:, :, i] * DqDx[:, :, j] / Qu)
                Fyy = torch.sum(DqDy[:, :, i] * DqDy[:, :, j] / Qu)

                # Fill symmetric entries simultaneously.
                F[ii,     jj    ] = Fxx;  F[jj,     ii    ] = Fxx
                F[ii,     jj + 1] = Fxy;  F[jj + 1, ii    ] = Fxy
                F[ii + 1, jj    ] = Fyx;  F[jj,     ii + 1] = Fyx
                F[ii + 1, jj + 1] = Fyy;  F[jj + 1, ii + 1] = Fyy

        return DtDxDy * I_mean * F

    # ================================================================
    # CASE KA — βm arbitrary, intensities known.  F = F11,  (2M × 2M).
    # Sun & Guan 2021, Eq (28) with arbitrary βi, βm.
    # ================================================================
    if ImType == 'KA':
        F = torch.zeros(2 * M, 2 * M, device=device, dtype=dtype)

        for i in range(M):
            for j in range(i, M):
                ii, jj = 2 * i, 2 * j
                bij = beta[i] * beta[j]

                Fxx = bij * torch.sum(DqDx[:, :, i] * DqDx[:, :, j] / Qu)
                Fxy = bij * torch.sum(DqDx[:, :, i] * DqDy[:, :, j] / Qu)
                Fyx = bij * torch.sum(DqDy[:, :, i] * DqDx[:, :, j] / Qu)
                Fyy = bij * torch.sum(DqDy[:, :, i] * DqDy[:, :, j] / Qu)

                F[ii,     jj    ] = Fxx;  F[jj,     ii    ] = Fxx
                F[ii,     jj + 1] = Fxy;  F[jj + 1, ii    ] = Fxy
                F[ii + 1, jj    ] = Fyx;  F[jj,     ii + 1] = Fyx
                F[ii + 1, jj + 1] = Fyy;  F[jj + 1, ii + 1] = Fyy

        return DtDxDy * I_mean * F

    # ================================================================
    # CASE UI — unknown common intensity I.  F = [[F11, F12],[F21, F22]],
    # size (2M+1, 2M+1).  UKI in the paper.
    # F11: Eq (28) with βm = 1.
    # F12 / F21: Eq (29) with βm = 1 → Eq (12) form.
    # F22: Eq (30) → Eq (13) form.
    # ================================================================
    if ImType == 'UI':
        F   = torch.zeros(2 * M + 1, 2 * M + 1, device=device, dtype=dtype)
        col = 2 * M   # index of the single unknown intensity I

        # F11 block — identical structure to KI.
        for i in range(M):
            for j in range(i, M):
                ii, jj = 2 * i, 2 * j

                Fxx = torch.sum(DqDx[:, :, i] * DqDx[:, :, j] / Qu)
                Fxy = torch.sum(DqDx[:, :, i] * DqDy[:, :, j] / Qu)
                Fyx = torch.sum(DqDy[:, :, i] * DqDx[:, :, j] / Qu)
                Fyy = torch.sum(DqDy[:, :, i] * DqDy[:, :, j] / Qu)

                F[ii,     jj    ] = Fxx;  F[jj,     ii    ] = Fxx
                F[ii,     jj + 1] = Fxy;  F[jj + 1, ii    ] = Fxy
                F[ii + 1, jj    ] = Fyx;  F[jj,     ii + 1] = Fyx
                F[ii + 1, jj + 1] = Fyy;  F[jj + 1, ii + 1] = Fyy

        F[:2 * M, :2 * M] *= DtDxDy * I_mean

        # F12 / F21 — Sun & Guan 2021, Eq (12) with Gaussian approx.
        # F(θij, I) = DtDxDy * Σk Q(k)/Qu(k) * ∂qi(k)/∂θij
        for i in range(M):
            ii = 2 * i
            Fx = DtDxDy * torch.sum(DqDx[:, :, i] * Q / Qu)
            Fy = DtDxDy * torch.sum(DqDy[:, :, i] * Q / Qu)
            F[ii,     col] = Fx;  F[col, ii    ] = Fx
            F[ii + 1, col] = Fy;  F[col, ii + 1] = Fy

        # F22 scalar — Sun & Guan 2021, Eq (13) with Gaussian approx.
        # F(I, I) = (DtDxDy / I_mean) * Σk Q^2(k) / Qu(k)
        F[col, col] = (DtDxDy / I_mean) * torch.sum(Q * Q / Qu)

        return F

    # ================================================================
    # CASE UA — all Im unknown and arbitrary.  F = [[F11, F12],[F21, F22]],
    # size (3M, 3M).  UAA in the paper.
    # Intensity index offset: Im occupies columns/rows 2M + m.
    # F11: Eq (28) with arbitrary βi, βm.
    # F12 / F21: Eq (29).
    # F22: Eq (30).
    # ================================================================
    if ImType == 'UA':
        F = torch.zeros(3 * M, 3 * M, device=device, dtype=dtype)

        # F11 block — identical structure to KA.
        for i in range(M):
            for j in range(i, M):
                ii, jj = 2 * i, 2 * j
                bij = beta[i] * beta[j]

                Fxx = bij * torch.sum(DqDx[:, :, i] * DqDx[:, :, j] / Qu)
                Fxy = bij * torch.sum(DqDx[:, :, i] * DqDy[:, :, j] / Qu)
                Fyx = bij * torch.sum(DqDy[:, :, i] * DqDx[:, :, j] / Qu)
                Fyy = bij * torch.sum(DqDy[:, :, i] * DqDy[:, :, j] / Qu)

                F[ii,     jj    ] = Fxx;  F[jj,     ii    ] = Fxx
                F[ii,     jj + 1] = Fxy;  F[jj + 1, ii    ] = Fxy
                F[ii + 1, jj    ] = Fyx;  F[jj,     ii + 1] = Fyx
                F[ii + 1, jj + 1] = Fyy;  F[jj + 1, ii + 1] = Fyy

        F[:2 * M, :2 * M] *= DtDxDy * I_mean

        # F12 / F21 — Sun & Guan 2021, Eq (29).
        # F(θij, Im) = DtDxDy * βi * Σk qm(k)/Qu(k) * ∂qi(k)/∂θij
        # Note: i indexes the location emitter, j indexes the intensity emitter.
        for i in range(M):
            ii = 2 * i
            bi = beta[i]
            for j in range(M):
                jj = 2 * M + j
                Fx = DtDxDy * bi * torch.sum(DqDx[:, :, i] * q_all[:, :, j] / Qu)
                Fy = DtDxDy * bi * torch.sum(DqDy[:, :, i] * q_all[:, :, j] / Qu)
                F[ii,     jj] = Fx;  F[jj, ii    ] = Fx
                F[ii + 1, jj] = Fy;  F[jj, ii + 1] = Fy

        # F22 — Sun & Guan 2021, Eq (30).
        # F(Ii, Im) = (DtDxDy / I_mean) * Σk qi(k) qm(k) / Qu(k)
        for i in range(M):
            for j in range(i, M):
                ii = 2 * M + i
                jj = 2 * M + j
                val = (DtDxDy / I_mean) * torch.sum(
                    q_all[:, :, i] * q_all[:, :, j] / Qu
                )
                F[ii, jj] = val
                F[jj, ii] = val

        return F

    raise ValueError(f"Invalid ImType '{ImType}'.")


def gauss2d_ugia_f_torch(
    system: SystemConfig,
    ImType: str,
):
    """UGIA-F estimator for 2D Gaussian PSF.

    Generates one unbiased Gaussian sample that achieves the Cramér–Rao
    lower bound (CRLB) on emitter *positions*, given the intensity case
    specified by ImType.  Used as a theoretical benchmark for practical
    localization algorithms.

    The sample is drawn as

        xyF = xy_true + W,   W ~ N(0, CRLB_xy)

    where CRLB_xy is the (2M × 2M) position covariance that a CRLB-
    achieving estimator would produce.  The FWHM CRLB for each coordinate
    follows from the diagonal: FWHM_CRLB_i = 2*sqrt(2*ln2) * sqrt(CRLB_xy[i,i]).

    Parameters
    ----------
    system : SystemConfig
        Must contain Gaussian2DPSF, 2-D EmitterData (frame 0 used),
        Camera, and Noise.
    ImType : str
        'KI', 'KA', 'UI', 'UA' — same convention as gauss2d_fisher_torch.

    Returns
    -------
    xyF      : (M, 2) tensor
        UGIA-F location estimates for this frame.
    F_xy     : (2M, 2M) tensor
        F11 block — position–position sub-matrix of the full FIM.
        For KI/KA this equals the full FIM; for UI/UA it is a sub-block
        and does *not* by itself give the CRLB for locations.
    crlb_xy  : (2M, 2M) tensor
        Position CRLB matrix — the (2M × 2M) upper-left block of the
        *full* FIM inverse.  For KI/KA: crlb_xy = pinv(F_xy).
        For UI/UA: crlb_xy = F_full_inv[:2M, :2M], which accounts for
        the Schur-complement degradation from unknown intensities and is
        strictly larger than pinv(F_xy).

    Notes
    -----
    Reference: Y. Sun & Y. Guan, JOSA A 38, 1830 (2021), and
    Y. Sun, Sci. Rep. 8, 17211 (2018) for the UGIA estimator definition.
    """

    if ImType not in ('KI', 'KA', 'UI', 'UA'):
        raise ValueError(
            f"ImType must be 'KI', 'KA', 'UI', or 'UA'; got '{ImType}'."
        )

    xy     = system.emitter.positions[0]   # (M, 2)
    device = xy.device
    dtype  = xy.dtype
    M      = xy.shape[0]

    if M < 1:
        raise ValueError("No emitters in frame.")

    # ------------------------------------------------------------------
    # Step 1: Full Fisher information matrix.
    # ------------------------------------------------------------------
    F_full = gauss2d_fisher_torch(system, ImType)

    # ------------------------------------------------------------------
    # Step 2: F11 — the (2M × 2M) position–position block.
    # Returned for reference; do NOT invert this for the CRLB when
    # intensities are unknown (UI / UA) — see Step 3.
    # ------------------------------------------------------------------
    F_xy = F_full[: 2 * M, : 2 * M].clone()

    # ------------------------------------------------------------------
    # Step 3: Position CRLB = upper-left (2M × 2M) block of F_full⁻¹.
    #
    # For KI / KA  : F_full IS F11,  so F_full_inv[:2M, :2M] = pinv(F_xy).
    # For UI / UA  : intensities are nuisance parameters.  The correct
    #                position CRLB is the Schur complement
    #                  (F11 - F12 F22⁻¹ F21)⁻¹
    #                obtained directly as F_full_inv[:2M, :2M].
    #                This is strictly ≥ pinv(F11) — the "intensity penalty"
    #                demonstrated in Sun & Guan (2021).
    # ------------------------------------------------------------------
    F_full_inv = torch.linalg.pinv(F_full)               # full inverse
    crlb_xy    = F_full_inv[: 2 * M, : 2 * M]            # (2M, 2M)

    # ------------------------------------------------------------------
    # Step 4: Draw W ~ N(0, crlb_xy).
    #
    # Eigendecomposition of the symmetric PSD matrix crlb_xy:
    #   crlb_xy = V @ diag(λ) @ V^T
    # Square-root factor: L = V @ diag(sqrt(λ))  →  L @ L^T = crlb_xy.
    # Eigenvalues are clamped to [0, ∞) for numerical safety; near-zero
    # values correspond to directions with very high CRLB (near-infinite
    # uncertainty) and contribute negligible noise.
    # ------------------------------------------------------------------
    eigenvalues, V = torch.linalg.eigh(crlb_xy)          # ascending order
    eigenvalues    = torch.clamp(eigenvalues, min=0.0)

    # L[:, j] = V[:, j] * sqrt(λ_j)  — broadcast scales each column.
    L = V * torch.sqrt(eigenvalues)                       # (2M, 2M)

    z = torch.randn(2 * M, device=device, dtype=dtype)    # (2M,)
    W = L @ z                                             # (2M,)  ~ N(0, crlb_xy)

    # ------------------------------------------------------------------
    # Step 5: UGIA-F sample.
    # ------------------------------------------------------------------
    xyF = xy + W.view(M, 2)                               # (M, 2)

    return xyF, F_xy, crlb_xy


def gauss2d_ugia_f_batch_torch(
    xy_all:  torch.Tensor,
    Im0:     float,
    psf,
    camera,
    noise,
    ImType:  str = 'KI',
) -> torch.Tensor:
    """Batched UGIA-F estimator — processes N frames simultaneously on GPU.

    Equivalent to calling gauss2d_ugia_f_torch N times but fully
    vectorised: all Fisher matrices are assembled and inverted in one
    batched GPU operation, giving ~10-50× speedup over the sequential loop.

    Supports ImType='KI' only (the common benchmark case).

    Parameters
    ----------
    xy_all : (N, M, 2) tensor — true emitter positions for all N frames
    Im0    : float            — emitter intensity (photons/s), same for all
    psf    : Gaussian2DPSF
    camera : Camera
    noise  : Noise
    ImType : str, must be 'KI'

    Returns
    -------
    xyF_all : (N, M, 2) tensor — UGIA-F estimates for all N frames
    """
    import math as _math

    if ImType != 'KI':
        raise NotImplementedError(
            "gauss2d_ugia_f_batch_torch supports ImType='KI' only."
        )

    N, M, _ = xy_all.shape
    device   = xy_all.device
    dtype    = xy_all.dtype

    sigma = psf.sigma
    Kx, Ky = camera.Kx, camera.Ky
    Dx, Dy = camera.Dx, camera.Dy
    Dt     = camera.Dt

    b = noise.b.to(device=device, dtype=dtype)   # (Ky, Kx)
    G = noise.G.to(device=device, dtype=dtype)   # (Ky, Kx)

    kx = torch.arange(Kx, device=device, dtype=dtype)
    ky = torch.arange(Ky, device=device, dtype=dtype)

    Dxk0 = Dx * kx          # (Kx,)
    Dxk1 = Dx * (kx + 1)
    Dyk0 = Dy * ky           # (Ky,)
    Dyk1 = Dy * (ky + 1)

    sqrt2pi_sigma = _math.sqrt(2.0 * _math.pi) * sigma
    I_mean        = float(Im0)
    DtDxDy        = Dt * Dx * Dy

    # ------------------------------------------------------------------
    # Per-emitter PSF and derivatives — batched over N frames and M emitters
    # xy_all: (N, M, 2)
    # x0, y0: (N, M)
    # ------------------------------------------------------------------
    x0 = xy_all[:, :, 0]   # (N, M)
    y0 = xy_all[:, :, 1]   # (N, M)

    # Pixel-edge offsets: (N, M, Kx) and (N, M, Ky)
    # Dxk0: (Kx,) → broadcast with x0: (N, M, 1)
    Dxk0_n = (Dxk0.view(1, 1, Kx) - x0.unsqueeze(2)) / sigma   # (N, M, Kx)
    Dxk1_n = (Dxk1.view(1, 1, Kx) - x0.unsqueeze(2)) / sigma
    Dyk0_n = (Dyk0.view(1, 1, Ky) - y0.unsqueeze(2)) / sigma   # (N, M, Ky)
    Dyk1_n = (Dyk1.view(1, 1, Ky) - y0.unsqueeze(2)) / sigma

    # PSF marginals: (N, M, Kx) and (N, M, Ky)
    qx    = (Qfunc_torch(-Dxk1_n) - Qfunc_torch(-Dxk0_n)) / Dx
    qy    = (Qfunc_torch(-Dyk1_n) - Qfunc_torch(-Dyk0_n)) / Dy
    DqxDx = (torch.exp(-Dxk0_n**2 / 2) - torch.exp(-Dxk1_n**2 / 2)) / (Dx * sqrt2pi_sigma)
    DqyDy = (torch.exp(-Dyk0_n**2 / 2) - torch.exp(-Dyk1_n**2 / 2)) / (Dy * sqrt2pi_sigma)

    # 2-D PSF and derivatives: (N, M, Ky, Kx)
    # qxy[n,m,ky,kx] = qy[n,m,ky] * qx[n,m,kx]
    qxy  = qy.unsqueeze(3) * qx.unsqueeze(2)           # (N, M, Ky, Kx)
    dqdx = qy.unsqueeze(3) * DqxDx.unsqueeze(2)        # (N, M, Ky, Kx)
    dqdy = DqyDy.unsqueeze(3) * qx.unsqueeze(2)        # (N, M, Ky, Kx)

    # Q = Σm qm(k) for KI (beta=1): (N, Ky, Kx)
    Q = qxy.sum(dim=1)                                  # (N, Ky, Kx)

    # Qu = Q + (b+G)/I_mean: (N, Ky, Kx)
    bg = ((b + G) / I_mean).unsqueeze(0)                # (1, Ky, Kx)
    Qu = torch.clamp(Q + bg, min=1e-10)                 # (N, Ky, Kx)

    # ------------------------------------------------------------------
    # Build Fisher matrix F: (N, 2M, 2M) for KI case
    # F[n, 2i, 2j]     = DtDxDy * I_mean * Σk DqDx[n,i,k] * DqDx[n,j,k] / Qu[n,k]
    # ------------------------------------------------------------------
    # dqdx/Qu: (N, M, Ky, Kx)
    dqdx_Qu = dqdx / Qu.unsqueeze(1)                    # (N, M, Ky, Kx)
    dqdy_Qu = dqdy / Qu.unsqueeze(1)                    # (N, M, Ky, Kx)

    # Stack into (N, 2M, Ky, Kx): interleave x and y derivatives
    # grad[n, 2i,   :, :] = dqdx[n, i, :, :]
    # grad[n, 2i+1, :, :] = dqdy[n, i, :, :]
    grad_num = torch.zeros(N, 2*M, Ky, Kx, device=device, dtype=dtype)
    grad_den = torch.zeros(N, 2*M, Ky, Kx, device=device, dtype=dtype)
    for i in range(M):
        grad_num[:, 2*i,   :, :] = dqdx[:, i, :, :]
        grad_num[:, 2*i+1, :, :] = dqdy[:, i, :, :]
        grad_den[:, 2*i,   :, :] = dqdx_Qu[:, i, :, :]
        grad_den[:, 2*i+1, :, :] = dqdy_Qu[:, i, :, :]

    # F[n, a, b] = DtDxDy * I_mean * Σ_{ky,kx} grad_num[n,a] * grad_den[n,b]
    # Reshape to (N, 2M, Ky*Kx) for batched matmul
    gn = grad_num.view(N, 2*M, Ky*Kx)   # (N, 2M, K)
    gd = grad_den.view(N, 2*M, Ky*Kx)   # (N, 2M, K)

    # F = DtDxDy * I_mean * gn @ gd^T  (batched)
    F = DtDxDy * I_mean * torch.bmm(gn, gd.transpose(1, 2))   # (N, 2M, 2M)

    # ------------------------------------------------------------------
    # Batched pinv and eigh — fully on GPU
    # ------------------------------------------------------------------
    F_inv      = torch.linalg.pinv(F)                    # (N, 2M, 2M)
    crlb_xy    = F_inv                                   # KI: full inv = crlb

    eigenvalues, V = torch.linalg.eigh(crlb_xy)         # (N, 2M), (N, 2M, 2M)
    eigenvalues    = torch.clamp(eigenvalues, min=0.0)

    # L[n, :, j] = V[n, :, j] * sqrt(λ[n,j])
    L = V * torch.sqrt(eigenvalues).unsqueeze(1)         # (N, 2M, 2M)

    # Draw W[n] ~ N(0, crlb_xy[n])
    z = torch.randn(N, 2*M, 1, device=device, dtype=dtype)   # (N, 2M, 1)
    W = torch.bmm(L, z).squeeze(2)                       # (N, 2M)

    # UGIA-F sample: xyF[n] = xy_all[n] + W[n].view(M, 2)
    xyF_all = xy_all + W.view(N, M, 2)                   # (N, M, 2)

    return xyF_all


"""
EM localization functions for the 2D Gaussian PSF.
Append to Gauss2D.py.

Reference: Sun & Fu (2014, QBI) for the EM algorithm;
           Sun (2013, JBO) for the imaging model and Gaussian-to-Poisson
           approximation (Sec. 2.2.3).
"""

# ============================================================
# Single-Emitter Maximum Likelihood — Newton method
# ============================================================

def gauss2d_seml_newton(
    system: SystemConfig,
    xy0: torch.Tensor,
    alpha: float,
    iter_no: int,
) -> tuple:
    """Single-Emitter Maximum Likelihood (SEML) via Newton's method.

    Maximises the per-pixel Poisson log-likelihood for one 2D Gaussian
    emitter.  The Gaussian-to-Poisson noise correction
        S(k) = U(k) + Dt*Dx*Dy*(G(k) - mu(k)) * mask_0(k)
    is applied internally (Sun 2013, JBO Sec. 2.2.3).

    When called from :func:`gauss2d_em_torch` the noise is zeroed
    (b=mu=G=0), so the correction evaluates to zero and S = Sm
    (the EM-attributed image) is used directly.

    Parameters
    ----------
    system : SystemConfig
        Must contain Gaussian2DPSF, Camera, per-pixel Noise, and
        ``frames`` set to the raw (or attributed) image, shape (1, Ky, Kx).
    xy0 : (2,) tensor
        Initial (x, y) in nm.  NaN triggers centroid initialisation.
    alpha : float
        Newton step scale (1.0 = standard; >1 can accelerate).
    iter_no : int
        Number of Newton iterations.

    Returns
    -------
    xy : (iter_no+1, 2) tensor — position at each iteration.
    I  : float               — emitter intensity used (photons/s).
    """
    psf    = system.psf
    camera = system.camera
    noise  = system.noise

    if psf.psf_type != "gaussian2d":
        raise ValueError("gauss2d_seml_newton requires Gaussian2DPSF.")
    if system.frames is None:
        raise ValueError("SystemConfig.frames must be set (shape 1, Ky, Kx).")

    # Device / dtype follow the initial position estimate.
    device = xy0.device
    dtype  = xy0.dtype

    sigma  = psf.sigma
    Kx, Ky = camera.Kx, camera.Ky
    Dx, Dy = camera.Dx, camera.Dy
    Dt     = camera.Dt

    # Align all tensors to the same device / dtype.
    # The frame is already Gaussian-corrected by the caller.
    U  = system.frames[0].to(device=device, dtype=dtype)   # (Ky, Kx)
    b  = noise.b.to(device=device, dtype=dtype)             # (Ky, Kx)
    G  = noise.G.to(device=device, dtype=dtype)             # (Ky, Kx)

    I0     = float(system.emitter.Im[0, 0].item())
    DtDxDy = Dt * Dx * Dy

    # ------------------------------------------------------------------
    # Gaussian-to-Poisson correction (Sun 2013, JBO Sec. 2.2.3):
    #   S(k) = U(k) + Dt*Dx*Dy*(G(k) - mu(k)) * mask_0(k)
    # Applied to non-zero pixels only (mask_0).  When called standalone
    # this converts the raw camera frame to an equivalent Poisson
    # observation.  When called from gauss2d_em_torch the noise is zeroed
    # (b=mu=G=0) so the correction evaluates to zero and S = U = Sm.
    # ------------------------------------------------------------------
    mu     = noise.mu.to(device=device, dtype=dtype)          # (Ky, Kx)
    mask_0 = (U != 0).to(dtype)                               # (Ky, Kx)
    S      = U + DtDxDy * (G - mu) * mask_0

    # Per-pixel background photon count  [photons / pixel / frame].
    s_bg = DtDxDy * (b + G)           # (Ky, Kx)

    # Estimate I0 from total photon count over non-zero pixels if needed.
    if math.isnan(I0):
        total_bg = (s_bg * mask_0.to(dtype)).sum().item()
        I0       = (S.sum().item() - total_bg) / DtDxDy

    I       = I0
    DtDxDyI = DtDxDy * I

    # Initialise position from brightest pixel if NaN.
    if torch.isnan(xy0).any():
        flat_idx = torch.argmax(U)
        kymax    = int(flat_idx // Kx)
        kxmax    = int(flat_idx %  Kx)
        xy0      = torch.tensor(
            [Dx * (kxmax + 0.5), Dy * (kymax + 0.5)],
            device=device, dtype=dtype,
        )

    # Allocate storage.
    xy    = torch.zeros(iter_no + 1, 2, device=device, dtype=dtype)
    Lxy   = torch.zeros(2,           device=device, dtype=dtype)
    Hxy   = torch.zeros(2, 2,        device=device, dtype=dtype)
    xy[0] = xy0

    kx = torch.arange(Kx, device=device, dtype=dtype)
    ky = torch.arange(Ky, device=device, dtype=dtype)
    sqrt2pi_sigma = math.sqrt(2.0 * math.pi) * sigma

    for j in range(iter_no):
        x0 = xy0[0];  y0 = xy0[1]

        Dxk1 = (Dx * (kx + 1.0) - x0) / sigma
        Dxk0 = (Dx *  kx        - x0) / sigma
        Dyk1 = (Dy * (ky + 1.0) - y0) / sigma
        Dyk0 = (Dy *  ky        - y0) / sigma

        # Pixel-averaged PSF marginals  [1/nm]  (Sun 2013, Eqs 45-46).
        qx = (Qfunc_torch(-Dxk1) - Qfunc_torch(-Dxk0)) / Dx   # (Kx,)
        qy = (Qfunc_torch(-Dyk1) - Qfunc_torch(-Dyk0)) / Dy   # (Ky,)

        # PSF x / y derivatives  [1/nm^3]  (Sun 2013, Eqs 37-38).
        DqxDx = (torch.exp(-0.5 * Dxk0**2) - torch.exp(-0.5 * Dxk1**2)) \
                / (Dx * sqrt2pi_sigma)                          # (Kx,)
        DqyDy = (torch.exp(-0.5 * Dyk0**2) - torch.exp(-0.5 * Dyk1**2)) \
                / (Dy * sqrt2pi_sigma)                          # (Ky,)

        q    = torch.outer(qy,    qx)     # (Ky, Kx)  qm(k)
        dqdx = torch.outer(qy,    DqxDx)  # (Ky, Kx)  ∂qm/∂x
        dqdy = torch.outer(DqyDy, qx)    # (Ky, Kx)  ∂qm/∂y

        # Predicted mean per pixel: background + signal.
        s = torch.clamp(s_bg + DtDxDyI * q, min=1e-10)        # (Ky, Kx)

        # Score: gradient of the Poisson log-likelihood w.r.t. (x, y).
        Ss1    = S / s - 1.0
        Lxy[0] = DtDxDyI * torch.sum(Ss1 * dqdx)
        Lxy[1] = DtDxDyI * torch.sum(Ss1 * dqdy)

        # 2×2 Fisher information matrix as Hessian approximation
        # (Fisher scoring; Sun & Guan 2021, Eq 28, single-emitter case).
        Qrp = torch.clamp(q + (b + G) / I, min=1e-10)         # (Ky, Kx)
        Hxy[0, 0] = torch.sum(dqdx * dqdx / Qrp)
        Hxy[0, 1] = torch.sum(dqdx * dqdy / Qrp)
        Hxy[1, 0] = Hxy[0, 1]
        Hxy[1, 1] = torch.sum(dqdy * dqdy / Qrp)
        Hxy_scaled = DtDxDyI * Hxy

        # Newton step  Δxy = H^{-1} L  (regularised against singularity).
        H_reg = Hxy_scaled + 1e-6 * torch.eye(2, device=device, dtype=dtype)
        Dxy   = torch.linalg.solve(H_reg, Lxy)

        xy0       = xy0 + alpha * Dxy
        xy[j + 1] = xy0

    return xy, I


# ============================================================
# Single-Emitter Maximum Likelihood — gradient ascent
# ============================================================

def gauss2d_seml_gradient(
    system: SystemConfig,
    xy0: torch.Tensor,
    alpha: float,
    iter_no: int,
) -> tuple:
    """Single-Emitter Maximum Likelihood (SEML) via gradient ascent.

    Same model as :func:`gauss2d_seml_newton` but uses the raw gradient
    (no curvature correction).  The Gaussian-to-Poisson correction is
    applied internally; when called from EM the noise is zeroed so
    S = Sm directly.

    Parameters
    ----------
    system : SystemConfig
        Must contain Gaussian2DPSF, Camera, per-pixel Noise, and
        ``frames`` set to the raw (or attributed) image, shape (1, Ky, Kx).
    xy0 : (2,) tensor
        Initial (x, y) in nm.  NaN triggers centroid initialisation.
    alpha : float
        Gradient step scale.
    iter_no : int
        Number of gradient ascent iterations.

    Returns
    -------
    xy : (iter_no+1, 2) tensor
    I  : float
    """
    psf    = system.psf
    camera = system.camera
    noise  = system.noise

    if psf.psf_type != "gaussian2d":
        raise ValueError("gauss2d_seml_gradient requires Gaussian2DPSF.")
    if system.frames is None:
        raise ValueError("SystemConfig.frames must be set (shape 1, Ky, Kx).")

    device = xy0.device
    dtype  = xy0.dtype

    sigma  = psf.sigma
    Kx, Ky = camera.Kx, camera.Ky
    Dx, Dy = camera.Dx, camera.Dy
    Dt     = camera.Dt

    # The frame is already Gaussian-corrected by the caller.
    U  = system.frames[0].to(device=device, dtype=dtype)
    b  = noise.b.to(device=device, dtype=dtype)
    G  = noise.G.to(device=device, dtype=dtype)

    I0     = float(system.emitter.Im[0, 0].item())
    DtDxDy = Dt * Dx * Dy

    # Gaussian-to-Poisson correction — same as gauss2d_seml_newton.
    # S(k) = U(k) + Dt*Dx*Dy*(G(k) - mu(k)) * mask_0(k)
    mu     = noise.mu.to(device=device, dtype=dtype)
    mask_0 = (U != 0).to(dtype)
    S      = U + DtDxDy * (G - mu) * mask_0
    s_bg   = DtDxDy * (b + G)

    if math.isnan(I0):
        total_bg = (s_bg * mask_0.to(dtype)).sum().item()
        I0       = (S.sum().item() - total_bg) / DtDxDy

    I       = I0
    DtDxDyI = DtDxDy * I

    if torch.isnan(xy0).any():
        flat_idx = torch.argmax(U)
        kymax    = int(flat_idx // Kx)
        kxmax    = int(flat_idx %  Kx)
        xy0      = torch.tensor(
            [Dx * (kxmax + 0.5), Dy * (kymax + 0.5)],
            device=device, dtype=dtype,
        )

    xy    = torch.zeros(iter_no + 1, 2, device=device, dtype=dtype)
    Lxy   = torch.zeros(2,           device=device, dtype=dtype)
    xy[0] = xy0

    kx = torch.arange(Kx, device=device, dtype=dtype)
    ky = torch.arange(Ky, device=device, dtype=dtype)
    sqrt2pi_sigma = math.sqrt(2.0 * math.pi) * sigma

    for j in range(iter_no):
        x0 = xy0[0];  y0 = xy0[1]

        Dxk1 = (Dx * (kx + 1.0) - x0) / sigma
        Dxk0 = (Dx *  kx        - x0) / sigma
        Dyk1 = (Dy * (ky + 1.0) - y0) / sigma
        Dyk0 = (Dy *  ky        - y0) / sigma

        qx = (Qfunc_torch(-Dxk1) - Qfunc_torch(-Dxk0)) / Dx
        qy = (Qfunc_torch(-Dyk1) - Qfunc_torch(-Dyk0)) / Dy

        DqxDx = (torch.exp(-0.5 * Dxk0**2) - torch.exp(-0.5 * Dxk1**2)) \
                / (Dx * sqrt2pi_sigma)
        DqyDy = (torch.exp(-0.5 * Dyk0**2) - torch.exp(-0.5 * Dyk1**2)) \
                / (Dy * sqrt2pi_sigma)

        q    = torch.outer(qy,    qx)
        dqdx = torch.outer(qy,    DqxDx)
        dqdy = torch.outer(DqyDy, qx)

        s = torch.clamp(s_bg + DtDxDyI * q, min=1e-10)

        Ss1    = S / s - 1.0
        Lxy[0] = DtDxDyI * torch.sum(Ss1 * dqdx)
        Lxy[1] = DtDxDyI * torch.sum(Ss1 * dqdy)

        xy0       = xy0 + alpha * Lxy
        xy[j + 1] = xy0

    return xy, I


# ============================================================
# EM algorithm for multi-emitter localization
# ============================================================

def gauss2d_em_torch(
    system: SystemConfig,
    xy_init: torch.Tensor,
    GN: str = 'G',
    max_iter: int = 1000,
) -> torch.Tensor:
    """EM algorithm for multi-emitter localization with 2D Gaussian PSF.

    Implements the Expectation-Maximisation algorithm of Sun & Fu (2014,
    QBI).  Gaussian readout noise is absorbed into an equivalent Poisson
    background term (Sun 2013, JBO, Sec. 2.2.3).

    The Gaussian-to-Poisson correction is applied **once** to the raw
    observed frame at the start of this function:
        U(k) ← U(k) + Dt*Dx*Dy*(G(k) − mu(k))
    All subsequent computations (E-step attribution and SEML M-step) use
    this corrected frame, so no further correction is needed inside SEML.

    E-step : attribute photon counts to individual emitters proportional
             to each emitter's predicted PSF contribution.
    M-step : re-localise each emitter on its attributed image via SEML.

    Convergence: each emitter must achieve T=20 consecutive iterations
    with total position change |Δx|+|Δy| < 0.05 nm, or max_iter is hit.

    Parameters
    ----------
    system : SystemConfig
        Must contain Gaussian2DPSF, Camera, per-pixel Noise (b, mu, G),
        and ``frames`` with shape (1, Ky, Kx) for the raw observed frame.
    xy_init : (M, 2) tensor
        Initial emitter positions in nm.
    GN : str
        'N' → Newton SEML;  any other value → gradient-ascent SEML.
    max_iter : int
        Hard upper limit on EM outer iterations (safety net).

    Returns
    -------
    xy : (M, 2) tensor of estimated emitter positions in nm.
    """
    emitter = system.emitter
    psf     = system.psf
    camera  = system.camera
    noise   = system.noise

    if system.frames is None:
        raise ValueError("SystemConfig.frames must be set (shape 1, Ky, Kx).")

    # EM operates on a single frame.
    U      = system.frames[0].clone()  # clone so we don't modify the original
    device = U.device
    dtype  = U.dtype

    sigma  = psf.sigma
    Kx, Ky = camera.Kx, camera.Ky
    Dx, Dy = camera.Dx, camera.Dy
    Dt     = camera.Dt

    Im = emitter.Im[0]    # (M,) — intensities for frame 0
    M  = xy_init.shape[0]

    if M == 0:
        return torch.zeros(0, 2, device=device, dtype=dtype)

    # Align noise and initial positions to working device / dtype.
    b        = noise.b.to(device=device, dtype=dtype)    # (Ky, Kx)
    G        = noise.G.to(device=device, dtype=dtype)    # (Ky, Kx)
    mu       = noise.mu.to(device=device, dtype=dtype)   # (Ky, Kx)
    xy_init  = xy_init.to(device=device, dtype=dtype)

    # Gaussian-to-Poisson correction applied ONCE to the full frame
    # before the E-step — matching MATLAB Gauss2D_EM.m:
    #   U = double(U + Dx*Dy*Dt*(G-mu))
    # Applied unconditionally to ALL pixels (no mask_0), unlike the
    # masked correction inside standalone SEML.
    DxDyDt = Dx * Dy * Dt
    U = U + DxDyDt * (G - mu)    # (Ky, Kx)

    # ------------------------------------------------------------------
    # EM hyper-parameters (matching MATLAB reference implementation).
    # ------------------------------------------------------------------
    alpha_newton   = 1.0    # step size for Newton method is 1 since it already uses curvature information
    iter_no_newton = 1      # one step is sufficient for SEML  
    alpha_gradient = 2.0    # This pair of alpha and iterations for Gradient ascent is optimum by testing 
    iter_no_grad   = 3

    err_tol = 5e-2   # nm  — convergence threshold (|Δx|+|Δy| per iter)
    T       = 20     # consecutive converged iterations to declare done

    kx = torch.arange(Kx, device=device, dtype=dtype)
    ky = torch.arange(Ky, device=device, dtype=dtype)

    errFlg = torch.zeros(M, device=device, dtype=torch.int32)
    xy     = torch.zeros_like(xy_init)
    xy0    = xy_init.clone()

    sqrt2pi_sigma = math.sqrt(2.0 * math.pi) * sigma
    # (M,1) / (1,Kx) broadcast shapes for vectorised PSF over M emitters
    kx_v = kx.unsqueeze(0)   # (1, Kx)
    ky_v = ky.unsqueeze(0)   # (1, Ky)
    Im_v = Im.view(M, 1)     # (M, 1) — per-emitter intensity

    iter_count = 0
    while (errFlg < T).any() and iter_count < max_iter:
        iter_count += 1

        # --------------------------------------------------------------
        # Vectorised PSF for all M emitters simultaneously.
        # xy0: (M,2) → x0:(M,1), y0:(M,1) for broadcasting over pixels.
        # Outputs: sm_all (M, Ky, Kx), q_all (M, Ky, Kx),
        #          dqdx_all (M, Ky, Kx), dqdy_all (M, Ky, Kx)
        # --------------------------------------------------------------
        x0v = xy0[:, 0].view(M, 1)   # (M, 1)
        y0v = xy0[:, 1].view(M, 1)   # (M, 1)

        # PSF marginals — shape (M, Kx) and (M, Ky)
        Dxk1 = (Dx * (kx_v + 1.0) - x0v) / sigma   # (M, Kx)
        Dxk0 = (Dx *  kx_v        - x0v) / sigma
        Dyk1 = (Dy * (ky_v + 1.0) - y0v) / sigma   # (M, Ky)
        Dyk0 = (Dy *  ky_v        - y0v) / sigma

        qx_all    = (Qfunc_torch(-Dxk1) - Qfunc_torch(-Dxk0)) / Dx  # (M, Kx)
        qy_all    = (Qfunc_torch(-Dyk1) - Qfunc_torch(-Dyk0)) / Dy  # (M, Ky)
        DqxDx_all = (torch.exp(-0.5*Dxk0**2) - torch.exp(-0.5*Dxk1**2)) \
                    / (Dx * sqrt2pi_sigma)                            # (M, Kx)
        DqyDy_all = (torch.exp(-0.5*Dyk0**2) - torch.exp(-0.5*Dyk1**2)) \
                    / (Dy * sqrt2pi_sigma)                            # (M, Ky)

        # Outer products over pixels: (M, Ky, Kx)
        q_all    = qy_all.unsqueeze(2)    * qx_all.unsqueeze(1)
        dqdx_all = qy_all.unsqueeze(2)    * DqxDx_all.unsqueeze(1)
        dqdy_all = DqyDy_all.unsqueeze(2) * qx_all.unsqueeze(1)

        # Expected signal per pixel per emitter (M, Ky, Kx)
        sm_all = DxDyDt * Im_v.unsqueeze(2) * q_all

        # Predicted total frame: sum over emitters + background (Ky, Kx)
        u = torch.clamp(sm_all.sum(0) + DxDyDt * (b + G), min=1e-10)

        # --------------------------------------------------------------
        # E-step: attributed images for all M emitters  (M, Ky, Kx)
        # --------------------------------------------------------------
        Sm_all = sm_all * U.unsqueeze(0) / u.unsqueeze(0)

        # --------------------------------------------------------------
        # M-step: vectorised SEML over all M emitters (b=G=0).
        # Only update emitters that have not yet converged.
        # --------------------------------------------------------------
        active = (errFlg < T)   # (M,) bool mask

        if active.any():
            DtDxDyI_v = DxDyDt * Im_v.squeeze(1)  # (M,)

            if GN == 'N':
                for _ in range(iter_no_newton):
                    s_all  = torch.clamp(DtDxDyI_v.view(M,1,1) * q_all, min=1e-10)
                    Ss1    = Sm_all / s_all - 1.0                   # (M, Ky, Kx)
                    # Gradient (M,)
                    Lx = DtDxDyI_v * (Ss1 * dqdx_all).sum(dim=(1,2))
                    Ly = DtDxDyI_v * (Ss1 * dqdy_all).sum(dim=(1,2))
                    # FIM diagonal + cross term (M,)
                    Qrp = torch.clamp(q_all, min=1e-10)
                    F00 = DtDxDyI_v * (dqdx_all**2       / Qrp).sum(dim=(1,2))
                    F01 = DtDxDyI_v * (dqdx_all*dqdy_all / Qrp).sum(dim=(1,2))
                    F11 = DtDxDyI_v * (dqdy_all**2       / Qrp).sum(dim=(1,2))
                    det = F00 * F11 - F01 * F01 + 1e-30
                    dx  = alpha_newton * ( F11 * Lx - F01 * Ly) / det
                    dy  = alpha_newton * (-F01 * Lx + F00 * Ly) / det
                    # Update positions and re-compute PSF for next sub-iter
                    xy[:, 0] = torch.where(active, xy0[:, 0] + dx, xy0[:, 0])
                    xy[:, 1] = torch.where(active, xy0[:, 1] + dy, xy0[:, 1])
                    if _ < iter_no_newton - 1:
                        # Recompute PSF at updated positions for next sub-iter
                        x0v = xy[:, 0].view(M, 1)
                        y0v = xy[:, 1].view(M, 1)
                        Dxk1 = (Dx * (kx_v + 1.0) - x0v) / sigma
                        Dxk0 = (Dx *  kx_v        - x0v) / sigma
                        Dyk1 = (Dy * (ky_v + 1.0) - y0v) / sigma
                        Dyk0 = (Dy *  ky_v        - y0v) / sigma
                        qx_all    = (Qfunc_torch(-Dxk1) - Qfunc_torch(-Dxk0)) / Dx
                        qy_all    = (Qfunc_torch(-Dyk1) - Qfunc_torch(-Dyk0)) / Dy
                        DqxDx_all = (torch.exp(-0.5*Dxk0**2) - torch.exp(-0.5*Dxk1**2)) \
                                    / (Dx * sqrt2pi_sigma)
                        DqyDy_all = (torch.exp(-0.5*Dyk0**2) - torch.exp(-0.5*Dyk1**2)) \
                                    / (Dy * sqrt2pi_sigma)
                        q_all    = qy_all.unsqueeze(2)    * qx_all.unsqueeze(1)
                        dqdx_all = qy_all.unsqueeze(2)    * DqxDx_all.unsqueeze(1)
                        dqdy_all = DqyDy_all.unsqueeze(2) * qx_all.unsqueeze(1)
            else:
                x0v = xy0[:, 0].view(M, 1)
                y0v = xy0[:, 1].view(M, 1)
                for _ in range(iter_no_grad):
                    Dxk1 = (Dx * (kx_v + 1.0) - x0v) / sigma
                    Dxk0 = (Dx *  kx_v        - x0v) / sigma
                    Dyk1 = (Dy * (ky_v + 1.0) - y0v) / sigma
                    Dyk0 = (Dy *  ky_v        - y0v) / sigma
                    qx_all    = (Qfunc_torch(-Dxk1) - Qfunc_torch(-Dxk0)) / Dx
                    qy_all    = (Qfunc_torch(-Dyk1) - Qfunc_torch(-Dyk0)) / Dy
                    DqxDx_all = (torch.exp(-0.5*Dxk0**2) - torch.exp(-0.5*Dxk1**2)) \
                                / (Dx * sqrt2pi_sigma)
                    DqyDy_all = (torch.exp(-0.5*Dyk0**2) - torch.exp(-0.5*Dyk1**2)) \
                                / (Dy * sqrt2pi_sigma)
                    q_all    = qy_all.unsqueeze(2)    * qx_all.unsqueeze(1)
                    dqdx_all = qy_all.unsqueeze(2)    * DqxDx_all.unsqueeze(1)
                    dqdy_all = DqyDy_all.unsqueeze(2) * qx_all.unsqueeze(1)
                    s_all  = torch.clamp(DtDxDyI_v.view(M,1,1) * q_all, min=1e-10)
                    Ss1    = Sm_all / s_all - 1.0
                    Lx = DtDxDyI_v * (Ss1 * dqdx_all).sum(dim=(1,2))
                    Ly = DtDxDyI_v * (Ss1 * dqdy_all).sum(dim=(1,2))
                    x0v = x0v + (alpha_gradient * Lx).view(M, 1)
                    y0v = y0v + (alpha_gradient * Ly).view(M, 1)
                xy[:, 0] = torch.where(active, x0v.squeeze(1), xy0[:, 0])
                xy[:, 1] = torch.where(active, y0v.squeeze(1), xy0[:, 1])

        # Convergence check for all emitters at once.
        err_xy = (xy - xy0).abs().sum(dim=1)   # (M,)
        errFlg = torch.where(err_xy < err_tol, errFlg + 1, torch.zeros_like(errFlg))

        xy0 = xy.clone()

    return xy


def gauss2d_em_batch_torch(
    frames:  torch.Tensor,
    xy_init: torch.Tensor,
    system:  SystemConfig,
    GN:      str = 'N',
    n_iter:  int = 10,
    dtype:   torch.dtype = torch.float32,
) -> torch.Tensor:
    """Batched EM localization over B frames with fixed number of iterations.

    Processes all B frames simultaneously by adding a leading batch
    dimension to every tensor.  The convergence check of
    gauss2d_em_torch is replaced by a fixed n_iter outer EM iterations,
    which eliminates the Python while-loop and enables full GPU
    parallelism across frames.

    From empirical testing n_iter=8 already gives good accuracy;
    the default n_iter=10 provides a small safety margin.

    Parameters
    ----------
    frames : (B, Ky, Kx) tensor
        Raw observed frames (cast to float32 internally).
        The Gaussian-to-Poisson correction is applied inside this function.
    xy_init : (B, M, 2) tensor
        Initial emitter positions in nm for each of the B frames.
    system : SystemConfig
        Provides PSF, Camera, Noise, and emitter intensities.
        system.frames is ignored -- pass frames directly.
    GN : str
        'N' -> Newton SEML (default);  else -> gradient-ascent SEML.
    n_iter : int
        Number of EM outer iterations.  Default 10.

    Returns
    -------
    xy : (B, M, 2) tensor of estimated emitter positions in nm.
    """
    psf    = system.psf
    camera = system.camera
    noise  = system.noise

    sigma   = psf.sigma
    Kx, Ky  = camera.Kx, camera.Ky
    Dx, Dy  = camera.Dx, camera.Dy
    Dt      = camera.Dt
    DxDyDt  = Dx * Dy * Dt

    device = frames.device
    # dtype from argument: float32 default; pass float64 for high photon
    # counts (large N), where float32 pixel sums lose sub-nm precision.

    B = frames.shape[0]
    M = xy_init.shape[1]

    # Align noise to device / dtype
    b  = noise.b.to(device=device, dtype=dtype)    # (Ky, Kx)
    G  = noise.G.to(device=device, dtype=dtype)
    mu = noise.mu.to(device=device, dtype=dtype)

    # Gaussian-to-Poisson correction -- all B frames at once
    U = frames.to(dtype=dtype) + DxDyDt * (G - mu)   # (B, Ky, Kx)

    # Per-emitter intensity -- (M,) from system.emitter
    Im = system.emitter.Im[0].to(device=device, dtype=dtype)   # (M,)

    # EM hyper-parameters
    alpha_newton   = 1.0;  iter_no_newton = 1
    alpha_gradient = 1.0;  iter_no_grad   = 3

    sqrt2pi_sigma = math.sqrt(2.0 * math.pi) * sigma

    # Pixel grids -- shape (1, 1, Kx/Ky) for broadcasting over (B, M, pixels)
    kx_v = torch.arange(Kx, device=device, dtype=dtype).view(1, 1, Kx)
    ky_v = torch.arange(Ky, device=device, dtype=dtype).view(1, 1, Ky)
    Im_v = Im.view(1, M, 1)    # (1, M, 1)

    xy0 = xy_init.to(device=device, dtype=dtype).clone()   # (B, M, 2)
    xy  = xy0.clone()

    for _ in range(n_iter):

        # PSF for all B frames and M emitters.
        # x0v, y0v : (B, M, 1) -- broadcast over pixel axis.
        x0v = xy0[:, :, 0].unsqueeze(2)   # (B, M, 1)
        y0v = xy0[:, :, 1].unsqueeze(2)

        Dxk1 = (Dx * (kx_v + 1.0) - x0v) / sigma   # (B, M, Kx)
        Dxk0 = (Dx *  kx_v        - x0v) / sigma
        Dyk1 = (Dy * (ky_v + 1.0) - y0v) / sigma   # (B, M, Ky)
        Dyk0 = (Dy *  ky_v        - y0v) / sigma

        qx_all    = (Qfunc_torch(-Dxk1) - Qfunc_torch(-Dxk0)) / Dx
        qy_all    = (Qfunc_torch(-Dyk1) - Qfunc_torch(-Dyk0)) / Dy
        DqxDx_all = (torch.exp(-0.5*Dxk0**2) - torch.exp(-0.5*Dxk1**2)) \
                    / (Dx * sqrt2pi_sigma)
        DqyDy_all = (torch.exp(-0.5*Dyk0**2) - torch.exp(-0.5*Dyk1**2)) \
                    / (Dy * sqrt2pi_sigma)

        # (B, M, Ky, Kx)
        q_all    = qy_all.unsqueeze(3)    * qx_all.unsqueeze(2)
        dqdx_all = qy_all.unsqueeze(3)    * DqxDx_all.unsqueeze(2)
        dqdy_all = DqyDy_all.unsqueeze(3) * qx_all.unsqueeze(2)

        # Expected signal per pixel per emitter -- (B, M, Ky, Kx)
        sm_all = DxDyDt * Im_v.unsqueeze(3) * q_all

        # Predicted total frame -- sum over M + background -- (B, Ky, Kx)
        u = torch.clamp(
            sm_all.sum(dim=1) + DxDyDt * (b + G).unsqueeze(0),
            min=1e-10,
        )

        # E-step: attributed images -- (B, M, Ky, Kx)
        Sm_all = sm_all * U.unsqueeze(1) / u.unsqueeze(1)

        # M-step: SEML with zero noise (b=G=0 passed implicitly)
        DtDxDyI_v = (DxDyDt * Im).view(1, M)   # (1, M)

        if GN == 'N':
            s_all = torch.clamp(DtDxDyI_v.view(1,M,1,1) * q_all, min=1e-10)
            Ss1   = Sm_all / s_all - 1.0
            Lx    = DtDxDyI_v * (Ss1 * dqdx_all).sum(dim=(2,3))   # (B,M)
            Ly    = DtDxDyI_v * (Ss1 * dqdy_all).sum(dim=(2,3))
            Qrp   = torch.clamp(q_all, min=1e-10)
            F00   = DtDxDyI_v * (dqdx_all**2         / Qrp).sum(dim=(2,3))
            F01   = DtDxDyI_v * (dqdx_all * dqdy_all / Qrp).sum(dim=(2,3))
            F11   = DtDxDyI_v * (dqdy_all**2         / Qrp).sum(dim=(2,3))
            det   = F00 * F11 - F01 * F01 + 1e-30
            xy[:, :, 0] = xy0[:, :, 0] + alpha_newton * ( F11*Lx - F01*Ly) / det
            xy[:, :, 1] = xy0[:, :, 1] + alpha_newton * (-F01*Lx + F00*Ly) / det
        else:
            x0v = xy0[:, :, 0].unsqueeze(2)   # (B, M, 1)
            y0v = xy0[:, :, 1].unsqueeze(2)
            for _g in range(iter_no_grad):
                s_all = torch.clamp(DtDxDyI_v.view(1,M,1,1) * q_all, min=1e-10)
                Ss1   = Sm_all / s_all - 1.0
                Lx    = DtDxDyI_v * (Ss1 * dqdx_all).sum(dim=(2,3))
                Ly    = DtDxDyI_v * (Ss1 * dqdy_all).sum(dim=(2,3))
                # Diagonal-Fisher preconditioning: the raw gradient grows
                # with the summed-frame intensity (~N) and diverges; dividing
                # by the diagonal Fisher makes the step scale as 1/N and stay
                # stable at any photon count.
                Qrp   = torch.clamp(q_all, min=1e-10)
                F00g  = DtDxDyI_v * (dqdx_all**2 / Qrp).sum(dim=(2,3))
                F11g  = DtDxDyI_v * (dqdy_all**2 / Qrp).sum(dim=(2,3))
                # Trust region: cap each step and keep emitters on the
                # frame, so a drifting emitter cannot reach the q->0 region
                # where F00g->0 and the step would explode (the gradient
                # divergence seen at high N on GPU).
                step_max = sigma
                sx = (alpha_gradient * Lx / (F00g + 1e-30)).clamp(-step_max, step_max)
                sy = (alpha_gradient * Ly / (F11g + 1e-30)).clamp(-step_max, step_max)
                x0v = (x0v + sx.unsqueeze(2)).clamp(0.0, Kx * Dx)
                y0v = (y0v + sy.unsqueeze(2)).clamp(0.0, Ky * Dy)
                if _g < iter_no_grad - 1:
                    Dxk1 = (Dx * (kx_v + 1.0) - x0v) / sigma
                    Dxk0 = (Dx *  kx_v        - x0v) / sigma
                    Dyk1 = (Dy * (ky_v + 1.0) - y0v) / sigma
                    Dyk0 = (Dy *  ky_v        - y0v) / sigma
                    qx_all    = (Qfunc_torch(-Dxk1) - Qfunc_torch(-Dxk0)) / Dx
                    qy_all    = (Qfunc_torch(-Dyk1) - Qfunc_torch(-Dyk0)) / Dy
                    DqxDx_all = (torch.exp(-0.5*Dxk0**2) - torch.exp(-0.5*Dxk1**2)) \
                                / (Dx * sqrt2pi_sigma)
                    DqyDy_all = (torch.exp(-0.5*Dyk0**2) - torch.exp(-0.5*Dyk1**2)) \
                                / (Dy * sqrt2pi_sigma)
                    q_all    = qy_all.unsqueeze(3)    * qx_all.unsqueeze(2)
                    dqdx_all = qy_all.unsqueeze(3)    * DqxDx_all.unsqueeze(2)
                    dqdy_all = DqyDy_all.unsqueeze(3) * qx_all.unsqueeze(2)
            xy[:, :, 0] = x0v.squeeze(2)
            xy[:, :, 1] = y0v.squeeze(2)

        xy0 = xy.clone()

    return xy   # (B, M, 2)




# ============================================================
# Information-sufficient SNR (per emitter, per frame)
# ============================================================

def gauss2d_snr_torch(
    system: SystemConfig,
    chunk_size: Optional[int] = None,
) -> torch.Tensor:
    """Per-emitter information-sufficient SNR for the 2D Gaussian PSF.

    Returns the signal-to-noise ratio in decibel for every emitter in every
    frame of a movie, ``snr_db`` of shape ``(N, M)``, parallel to
    ``emitter.Im``. The definition is the 80%-power-effective SNR of
    Sun (2020), extended to unequal emitter intensities and per-pixel
    non-uniform noise:

        n_bar_m = sum_k q_m(k) (b_k + G_k) / sum_k q_m(k)            (Eq. 27)
        gamma_m = I_m / n_bar_m
        SNR_m   = 10 log10(gamma_m) - 20 log10(sigma) - 11.02         (Eq. 29)

    Here q_m(k) is the pixel-integrated 2D Gaussian PSF, the fraction of
    emitter m's photons in pixel k; b_k and G_k are the Poisson-background
    and Gaussian-readout variance densities in photons/(s*nm^2); I_m is the
    emitter intensity in photons/s; sigma is the PSF width in nm. The
    -11.02 dB constant equals 10 log10(eta_s) with the unit-variance Gaussian
    figure of merit eta_s = -rho / [2 pi ln(1 - rho)] = 0.0791 at the
    information-sufficient level rho = 0.80; the precise value is used here.

    q_m(k) is computed independently of ``gauss2d_frame_torch`` using the
    same erf-based pixel integration, so no existing routine is invoked.
    Inactive emitters (Im == 0) return -inf.

    Parameters
    ----------
    system : SystemConfig
        Must contain a Gaussian2DPSF, a 2-D EmitterData, a Camera, and a
        Noise. Tensors may live on any device; all are aligned to the
        emitter device internally.
    chunk_size : int, optional
        Process the N frames ``chunk_size`` at a time to bound the
        (Nc, M, Ky, Kx) intermediate. ``None`` (default) does all N at once.

    Returns
    -------
    snr_db : (N, M) tensor
        Per-emitter SNR in dB. For the per-frame range and mean over the
        active emitters use ``snr_db.amin(1)``, ``snr_db.amax(1)``,
        ``snr_db.mean(1)``.

    References
    ----------
    Sun, Y. "Information sufficient segmentation and signal-to-noise ratio
    in stochastic optical localization nanoscopy." Opt. Lett. 45(21),
    6102-6105 (2020).
    M. Sun and Y. Sun, Electron. Lett. 58(2), 58-60 (2022), Eq. (11).
    """
    emitter = system.emitter
    psf = system.psf
    camera = system.camera
    noise = system.noise

    if psf.psf_type != "gaussian2d":
        raise ValueError(
            f"gauss2d_snr_torch requires Gaussian2DPSF; got '{psf.psf_type}'."
        )

    xy = emitter.positions          # (N, M, 2)
    Im = emitter.Im                 # (N, M)
    N, M, _ = xy.shape

    device = xy.device
    dtype = xy.dtype

    Kx, Ky = camera.Kx, camera.Ky
    Dx, Dy = camera.Dx, camera.Dy
    sigma = float(psf.sigma)

    Im = Im.to(device=device, dtype=dtype)

    # Per-pixel noise variance density b + G in photons/(s*nm^2). These are
    # densities, not per-pixel counts: gamma = I / (b + G).
    n_density = (noise.b + noise.G).to(device=device, dtype=dtype)   # (Ky, Kx)

    sqrt2 = math.sqrt(2.0)

    # Pixel-edge grid (nm).
    kx = torch.arange(Kx, device=device, dtype=dtype).view(1, 1, 1, Kx)
    ky = torch.arange(Ky, device=device, dtype=dtype).view(1, 1, Ky, 1)
    Dxk0, Dxk1 = Dx * kx, Dx * (kx + 1)
    Dyk0, Dyk1 = Dy * ky, Dy * (ky + 1)

    # Figure of merit (Eqs. 28-29).
    rho = 0.80
    eta_s = -rho / (2.0 * math.pi * math.log(1.0 - rho))   # 0.07911...
    eta = eta_s / (sigma * sigma)

    snr_db = torch.empty((N, M), device=device, dtype=dtype)
    chunk = N if chunk_size is None else int(chunk_size)
    if chunk <= 0:
        raise ValueError("chunk_size must be a positive integer or None.")

    for n0 in range(0, N, chunk):
        n1 = min(n0 + chunk, N)
        Nc = n1 - n0

        x0 = xy[n0:n1, :, 0].view(Nc, M, 1, 1)
        y0 = xy[n0:n1, :, 1].view(Nc, M, 1, 1)

        # Pixel-integrated PSF q_m(k) = Phi(upper) - Phi(lower),
        # Phi the standard normal CDF via erf, computed independently.
        qx = 0.5 * (torch.erf((Dxk1 - x0) / (sigma * sqrt2))
                    - torch.erf((Dxk0 - x0) / (sigma * sqrt2)))   # (Nc, M, 1, Kx)
        qy = 0.5 * (torch.erf((Dyk1 - y0) / (sigma * sqrt2))
                    - torch.erf((Dyk0 - y0) / (sigma * sqrt2)))   # (Nc, M, Ky, 1)
        q = qy * qx                                                # (Nc, M, Ky, Kx)

        # PSF-weighted local noise density n_bar_m (Eq. 27).
        num = (q * n_density).sum(dim=(-2, -1))    # sum_k q_m(k) (b_k + G_k)
        den = q.sum(dim=(-2, -1))                  # sum_k q_m(k)
        n_bar = num / den                          # (Nc, M)

        gamma = Im[n0:n1] / n_bar                  # (Nc, M), nm^2
        snr_db[n0:n1] = 10.0 * torch.log10(eta * gamma)

    active = Im > 0
    return torch.where(active, snr_db, torch.full_like(snr_db, float("-inf")))