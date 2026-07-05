"""
Astig3D.py
Functions related to the 3D astigmatic Gaussian PSF.

Implements the data-frame model of Sun (2013, JBO) and Sun & Guan (2021,
JOSA A) for the 3D astigmatic Gaussian PSF (cylindrical-lens 3D nanoscopy),
with the same per-pixel-noise and per-frame variable emitter-locations
extensions adopted elsewhere in SMLM_Lib. This is the 3D analog of
``gauss2d_frame_torch`` in Gauss2D.py and the equivalent of the MATLAB
``AS3D_Frame`` used in Sun (2024, Optics Express), Sec. 4.1.

The astigmatic PSF differs from the 2D Gaussian only in that the lateral
widths depend on the emitter axial position z:

    sigma_x(z) = sigmax0 * sqrt(1 + (z+c)^2/d^2 + Ax(z+c)^3/d^3 + Bx(z+c)^4/d^4)
    sigma_y(z) = sigmay0 * sqrt(1 + (z-c)^2/d^2 + Ay(z-c)^3/d^3 + By(z-c)^4/d^4)

so the pixel integration is the same separable normal-CDF computation, but
with a per-emitter width that differs between x and y.
"""

import math
from typing import Optional

import torch

from .config import (
    EmitterData,
    Astigmatic3DPSF,
    Camera,
    Noise,
    SystemConfig,
)
from .Gauss2D import Qfunc_torch


# ============================================================
# Helpers
# ============================================================

def astig_sigma_xy(z: torch.Tensor, psf: Astigmatic3DPSF):
    """Per-emitter astigmatic widths sigma_x(z), sigma_y(z).

    Parameters
    ----------
    z   : tensor of axial positions in nm (any shape).
    psf : Astigmatic3DPSF.

    Returns
    -------
    (sigma_x, sigma_y) : tensors matching ``z`` in shape, dtype, device.

    Notes
    -----
    The radicands are clamped to a small positive floor for numerical
    safety; within the usual axial range they are positive and the clamp
    has no effect.
    """
    c, d = psf.c, psf.d
    up = z + c
    un = z - c
    rad_x = (1.0 + up**2 / d**2 + psf.Ax * up**3 / d**3 + psf.Bx * up**4 / d**4)
    rad_y = (1.0 + un**2 / d**2 + psf.Ay * un**3 / d**3 + psf.By * un**4 / d**4)
    sigma_x = psf.sigmax0 * torch.sqrt(rad_x.clamp_min(1e-12))
    sigma_y = psf.sigmay0 * torch.sqrt(rad_y.clamp_min(1e-12))
    return sigma_x, sigma_y


# ============================================================
# Frame simulator
# ============================================================

def astig3d_frame_torch(
    system: SystemConfig,
    chunk_size: Optional[int] = None,
) -> torch.Tensor:
    """Generate SMLM frames using the 3D astigmatic Gaussian PSF.

    Implements the per-pixel mean photon count

        v(n, ky, kx) = Dt * Dx * Dy * sum_m  Im[n, m] * q_m(n, ky, kx)
                       + Dt * Dx * Dy * b(ky, kx)

    where q_m is the pixel-averaged astigmatic PSF, evaluated via the
    standard normal CDF along x and y separately with per-emitter widths
    sigma_x(z_m) and sigma_y(z_m). Poisson photon counts are sampled from
    v, and Gaussian readout noise (per-pixel mean mu and variance G) is
    added. This mirrors ``gauss2d_frame_torch`` exactly aside from the
    z-dependent, axis-dependent widths.

    Parameters
    ----------
    system : SystemConfig
        Must contain an Astigmatic3DPSF, a 3-D EmitterData (xyz), a Camera,
        and a Noise. Emitter positions, intensities, and noise tensors may
        live on any devices; they are aligned to ``emitter.xyz``'s device.
    chunk_size : int, optional
        If given, frames are processed ``chunk_size`` at a time along the N
        axis to bound the (Nc, M, Ky*Kx) intermediate tensor. ``None``
        (default) processes all N frames at once.

    Returns
    -------
    U : (N, Ky, Kx) tensor
        Simulated frames in photoelectron units (float). Cast to uint16
        downstream for camera-style output.
    """
    emitter = system.emitter
    psf = system.psf
    camera = system.camera
    noise = system.noise

    # ------------------------------------------------------------
    # Validate PSF type
    # ------------------------------------------------------------
    if psf.psf_type != "astigmatic3d":
        raise ValueError(
            f"astig3d_frame_torch requires Astigmatic3DPSF; got '{psf.psf_type}'."
        )

    # ------------------------------------------------------------
    # Unpack
    # ------------------------------------------------------------
    xyz = emitter.positions         # (N, M, 3)
    Im = emitter.Im                 # (N, M)
    N, M, _ = xyz.shape

    device = xyz.device
    dtype = xyz.dtype

    Kx, Ky = camera.Kx, camera.Ky
    Dx, Dy = camera.Dx, camera.Dy
    Dt = camera.Dt

    # ------------------------------------------------------------
    # Per-pixel noise quantities — align to the emitter's device.
    # ------------------------------------------------------------
    b_pix  = noise.poisson_mean_per_pixel(camera).to(device=device, dtype=dtype)
    g_mean = noise.gaussian_mean_per_pixel(camera).to(device=device, dtype=dtype)
    g_var  = noise.gaussian_var_per_pixel(camera).to(device=device, dtype=dtype)
    g_std  = torch.sqrt(g_var)

    # ------------------------------------------------------------
    # Pixel grid boundaries (independent of frames/emitters)
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
        x0 = xyz[n0:n1, :, 0].view(Nc, M, 1, 1)
        y0 = xyz[n0:n1, :, 1].view(Nc, M, 1, 1)
        z0 = xyz[n0:n1, :, 2].view(Nc, M, 1, 1)

        # Per-emitter astigmatic widths (depend on z, differ in x and y).
        sx, sy = astig_sigma_xy(z0, psf)                   # each (Nc, M, 1, 1)

        # Normalised pixel-edge offsets (per-emitter, per-axis width).
        Dxk0_n = (Dxk0 - x0) / sx
        Dxk1_n = (Dxk1 - x0) / sx
        Dyk0_n = (Dyk0 - y0) / sy
        Dyk1_n = (Dyk1 - y0) / sy

        # Pixel-integrated PSF: Phi(upper) - Phi(lower) via Q(-x) = Phi(x).
        qx = Qfunc_torch(-Dxk1_n) - Qfunc_torch(-Dxk0_n)   # (Nc, M, 1, Kx)
        qy = Qfunc_torch(-Dyk1_n) - Qfunc_torch(-Dyk0_n)   # (Nc, M, Ky, 1)
        q  = qy * qx                                         # (Nc, M, Ky, Kx)

        # Weighted sum over emitters via bmm.
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

def _astig_sigma_and_dsigma2_dz(z, psf):
    """Per-emitter sigma_x, sigma_y and d(sigma^2)/dz for the astigmatic PSF.

    Returns (sx, sy, dsx2, dsy2), each matching ``z`` in shape. With
    zx = (z+c)/d, zy = (z-c)/d:
        sx^2 = sx0^2 (1 + zx^2 + Ax zx^3 + Bx zx^4)
        d(sx^2)/dz = sx0^2 (2 zx + 3 Ax zx^2 + 4 Bx zx^3) / d
    (and analogously for y with zy, Ay, By).
    """
    c, d = psf.c, psf.d
    zx = (z + c) / d
    zy = (z - c) / d
    sx = psf.sigmax0 * torch.sqrt((1.0 + zx**2 + psf.Ax*zx**3 + psf.Bx*zx**4).clamp_min(1e-12))
    sy = psf.sigmay0 * torch.sqrt((1.0 + zy**2 + psf.Ay*zy**3 + psf.By*zy**4).clamp_min(1e-12))
    dsx2 = psf.sigmax0**2 * (2*zx + 3*psf.Ax*zx**2 + 4*psf.Bx*zx**3) / d
    dsy2 = psf.sigmay0**2 * (2*zy + 3*psf.Ay*zy**2 + 4*psf.By*zy**3) / d
    return sx, sy, dsx2, dsy2


def astig3d_fisher_torch(system: SystemConfig, ImType: str) -> torch.Tensor:
    """Fisher Information Matrix (FIM) for the 3D astigmatic Gaussian PSF.

    Port of ``AS3D_Fisher_Im.m`` (Sun, 01/2021); the 3D analog of
    ``gauss2d_fisher_torch``. Implements the data-frame CRLB of Sun (2013,
    JBO) for the astigmatic PSF, with the unknown-intensity cases of
    Sun & Guan (2021, JOSA A). Gaussian noise is approximated as Poisson
    (its mean does not enter the FIM).

    The astigmatic widths depend on z, so each emitter contributes a
    location triple (x, y, z); the z-information comes from the differential
    response of sigma_x(z) and sigma_y(z) through

        dq/dz = (dqx/dz) qy + qx (dqy/dz),
        dqx/dz = [d(sx^2)/dz / (2 Dx sx^2 sqrt(2pi))]
                 * (Dxk0 e^{-Dxk0^2/2} - Dxk1 e^{-Dxk1^2/2}),

    the closed form of the (x^2-1)e^{-x^2/2} pixel integral in the MATLAB.

    Parameter ordering within F:
        locations : (x1, y1, z1, ..., xM, yM, zM)   — indices 0 .. 3M-1
        intensity : I            for UI  — index 3M
                    I1, .., IM   for UA  — indices 3M .. 4M-1

    Parameters
    ----------
    system : SystemConfig
        Must contain an Astigmatic3DPSF, a 3-D EmitterData (frame 0 is used —
        CRLB is a per-frame quantity), a Camera, and a Noise with per-pixel
        b and G.
    ImType : {'KI', 'KA', 'UI', 'UA'}
        KI/KA -> (3M, 3M); UI -> (3M+1, 3M+1); UA -> (4M, 4M). Same meaning
        as in gauss2d_fisher_torch.

    Returns
    -------
    F : tensor — the Fisher information matrix.
    """
    if ImType not in ('KI', 'KA', 'UI', 'UA'):
        raise ValueError("ImType must be 'KI', 'KA', 'UI', or 'UA'.")

    emitter, psf, camera, noise = (system.emitter, system.psf,
                                   system.camera, system.noise)
    if psf.psf_type != "astigmatic3d":
        raise ValueError("astig3d_fisher_torch requires Astigmatic3DPSF.")

    xyz = emitter.positions[0]      # (M, 3) nm  (frame 0)
    Im  = emitter.Im[0]             # (M,)  photons/s
    M   = xyz.shape[0]
    device, dtype = xyz.device, xyz.dtype

    Kx, Ky = camera.Kx, camera.Ky
    Dx, Dy = camera.Dx, camera.Dy
    Dt = camera.Dt
    DtDxDy = Dt * Dx * Dy
    sqrt2pi = math.sqrt(2.0 * math.pi)

    b = noise.b.to(device=device, dtype=dtype)      # (Ky, Kx)
    G = noise.G.to(device=device, dtype=dtype)      # (Ky, Kx)

    I_mean = Im.mean()
    beta = (Im / I_mean) if ImType in ('KA', 'UA') \
           else torch.ones(M, device=device, dtype=dtype)        # (M,)

    # ----- Per-emitter marginals & derivatives (vectorised over M) -----
    x0 = xyz[:, 0:1]; y0 = xyz[:, 1:2]; z = xyz[:, 2]            # (M,1),(M,1),(M,)
    sx, sy, dsx2, dsy2 = _astig_sigma_and_dsigma2_dz(z, psf)     # (M,)
    sx = sx.unsqueeze(1); sy = sy.unsqueeze(1)                   # (M,1)
    dsx2 = dsx2.unsqueeze(1); dsy2 = dsy2.unsqueeze(1)

    kx = torch.arange(Kx, device=device, dtype=dtype)
    ky = torch.arange(Ky, device=device, dtype=dtype)
    Dxk0 = (Dx*kx - x0) / sx;  Dxk1 = (Dx*(kx+1) - x0) / sx      # (M, Kx)
    Dyk0 = (Dy*ky - y0) / sy;  Dyk1 = (Dy*(ky+1) - y0) / sy      # (M, Ky)

    qx = (Qfunc_torch(-Dxk1) - Qfunc_torch(-Dxk0)) / Dx          # (M, Kx)
    qy = (Qfunc_torch(-Dyk1) - Qfunc_torch(-Dyk0)) / Dy          # (M, Ky)

    DqxDx = (torch.exp(-Dxk0**2/2) - torch.exp(-Dxk1**2/2)) / (Dx*sqrt2pi*sx)
    DqyDy = (torch.exp(-Dyk0**2/2) - torch.exp(-Dyk1**2/2)) / (Dy*sqrt2pi*sy)

    # closed form of integral_{Dk0}^{Dk1} (x^2-1) e^{-x^2/2} dx = a e^{-a^2/2} - b e^{-b^2/2}
    iF0x = Dxk0*torch.exp(-Dxk0**2/2) - Dxk1*torch.exp(-Dxk1**2/2)
    iF0y = Dyk0*torch.exp(-Dyk0**2/2) - Dyk1*torch.exp(-Dyk1**2/2)
    DqxDz = (dsx2 / (2*Dx*sx**2*sqrt2pi)) * iF0x                 # (M, Kx)
    DqyDz = (dsy2 / (2*Dy*sy**2*sqrt2pi)) * iF0y                 # (M, Ky)

    # 2-D PSF and derivatives via separability; (M, Ky, Kx), rows=y cols=x.
    q  = qy.unsqueeze(2) * qx.unsqueeze(1)
    gx = qy.unsqueeze(2) * DqxDx.unsqueeze(1)
    gy = DqyDy.unsqueeze(2) * qx.unsqueeze(1)
    gz = DqxDz.unsqueeze(1) * qy.unsqueeze(2) + qx.unsqueeze(1) * DqyDz.unsqueeze(2)

    # Effective denominator Qu = Q + (b+G)/I_mean  (per-pixel), Q = sum_m beta_m q_m.
    Q  = torch.einsum('m,mij->ij', beta, q)                     # (Ky, Kx)
    Qu = torch.clamp(Q + (b + G) / I_mean, min=1e-10)           # (Ky, Kx)

    P  = Ky * Kx
    Wf = (1.0 / Qu).reshape(P)                                  # (P,)
    gxf = gx.reshape(M, P); gyf = gy.reshape(M, P); gzf = gz.reshape(M, P)
    qf  = q.reshape(M, P);  Qf  = Q.reshape(P)

    # ----- Location block (3M x 3M): F_ab[i,j] = sum_p g_a_i g_b_j / Qu -----
    def blk(A, B):                       # (M, M)
        return (A * Wf) @ B.t()
    Sxx, Sxy, Sxz = blk(gxf, gxf), blk(gxf, gyf), blk(gxf, gzf)
    Syy, Syz, Szz = blk(gyf, gyf), blk(gyf, gzf), blk(gzf, gzf)
    Syx, Szx, Szy = Sxy.t(), Sxz.t(), Syz.t()

    bb = beta.unsqueeze(1) * beta.unsqueeze(0)                  # (M, M); =1 if KI/UI
    T = torch.empty(M, 3, M, 3, device=device, dtype=dtype)
    T[:, 0, :, 0] = Sxx*bb; T[:, 0, :, 1] = Sxy*bb; T[:, 0, :, 2] = Sxz*bb
    T[:, 1, :, 0] = Syx*bb; T[:, 1, :, 1] = Syy*bb; T[:, 1, :, 2] = Syz*bb
    T[:, 2, :, 0] = Szx*bb; T[:, 2, :, 1] = Szy*bb; T[:, 2, :, 2] = Szz*bb
    F11 = T.reshape(3*M, 3*M) * (DtDxDy * I_mean)

    if ImType in ('KI', 'KA'):
        return F11

    if ImType == 'UI':
        F = torch.zeros(3*M + 1, 3*M + 1, device=device, dtype=dtype)
        F[:3*M, :3*M] = F11
        col = 3*M
        QW = Qf * Wf                                            # (P,)
        Fx = DtDxDy * (gxf @ QW); Fy = DtDxDy * (gyf @ QW); Fz = DtDxDy * (gzf @ QW)
        loc = torch.stack([Fx, Fy, Fz], dim=1).reshape(3*M)     # interleave x,y,z
        F[:3*M, col] = loc; F[col, :3*M] = loc
        F[col, col] = (DtDxDy / I_mean) * torch.sum(Qf * Qf * Wf)
        return F

    # ImType == 'UA'
    F = torch.zeros(4*M, 4*M, device=device, dtype=dtype)
    F[:3*M, :3*M] = F11
    Cx = (gxf * Wf) @ qf.t(); Cy = (gyf * Wf) @ qf.t(); Cz = (gzf * Wf) @ qf.t()  # (M,M)
    bi = beta.unsqueeze(1)                                      # (M,1) row scaling
    F12 = torch.stack([Cx*bi, Cy*bi, Cz*bi], dim=1).reshape(3*M, M) * DtDxDy
    F[:3*M, 3*M:] = F12; F[3*M:, :3*M] = F12.t()
    F22 = ((qf * Wf) @ qf.t()) * (DtDxDy / I_mean)              # (M, M)
    F[3*M:, 3*M:] = F22
    return F


def astig3d_ugia_f_torch(system: SystemConfig, ImType: str):
    """UGIA-F estimator for the 3D astigmatic Gaussian PSF.

    Port of ``AS3D_UGIA_F.m`` (Sun, 02/2020), generalised to the four
    intensity cases and made the 3D analog of ``gauss2d_ugia_f_torch``.
    Draws one unbiased Gaussian sample of the emitter *locations* that
    achieves the Cramer-Rao lower bound (CRLB):

        xyzF = xyz_true + W,    W ~ N(0, CRLB_xyz)

    where CRLB_xyz is the (3M x 3M) position covariance a CRLB-achieving
    estimator would produce. Used as a theoretical localization benchmark.
    The per-coordinate CRLB follows from the diagonal, e.g.
    FWHM_CRLB_i = 2*sqrt(2*ln2) * sqrt(CRLB_xyz[i, i]).

    The MATLAB ``AS3D_UGIA_F.m`` (known identical intensity, locations only)
    corresponds to ``ImType='KI'``; it computes F = AS3D_Fisher, F_ = pinv(F),
    and draws W with covariance F_ via an SVD factor. For a symmetric PSD
    covariance, the eigendecomposition factor used here is equivalent.

    Parameters
    ----------
    system : SystemConfig
        Must contain Astigmatic3DPSF, 3-D EmitterData (frame 0 used),
        Camera, and Noise.
    ImType : {'KI', 'KA', 'UI', 'UA'}
        Same convention as astig3d_fisher_torch.

    Returns
    -------
    xyzF     : (M, 3) tensor — UGIA-F location estimates for this frame.
    F_xyz    : (3M, 3M) tensor — position-position block (F11) of the full
        FIM. For KI/KA this is the full FIM; for UI/UA it is a sub-block and
        must NOT be inverted on its own for the CRLB (see crlb_xyz).
    crlb_xyz : (3M, 3M) tensor — position CRLB, the upper-left (3M x 3M)
        block of pinv(F_full). For KI/KA this equals pinv(F_xyz); for UI/UA
        it is the Schur complement (F11 - F12 F22^{-1} F21)^{-1}, strictly
        larger than pinv(F11) — the intensity penalty of Sun & Guan (2021).

    Notes
    -----
    References: Sun (2013, JBO) for the data-frame CRLB; Sun (2018, Sci.
    Rep.) for the UGIA estimator; Sun & Guan (2021, JOSA A) for the
    unknown-intensity cases.
    """
    if ImType not in ('KI', 'KA', 'UI', 'UA'):
        raise ValueError("ImType must be 'KI', 'KA', 'UI', or 'UA'.")

    xyz = system.emitter.positions[0]        # (M, 3)
    device, dtype = xyz.device, xyz.dtype
    M = xyz.shape[0]
    if M < 1:
        raise ValueError("No emitters in frame.")

    # Step 1: full Fisher information matrix for the requested intensity case.
    F_full = astig3d_fisher_torch(system, ImType)

    # Step 2: (3M x 3M) position-position block (returned for reference).
    F_xyz = F_full[:3 * M, :3 * M].clone()

    # Step 3: position CRLB = upper-left (3M x 3M) block of the full inverse.
    #   KI/KA: F_full is F11, so this equals pinv(F_xyz).
    #   UI/UA: Schur complement w.r.t. the unknown intensities (>= pinv(F11)).
    crlb_xyz = torch.linalg.pinv(F_full)[:3 * M, :3 * M]

    # Step 4: draw W ~ N(0, crlb_xyz) via the symmetric PSD square-root factor.
    eigenvalues, V = torch.linalg.eigh(crlb_xyz)          # ascending
    eigenvalues = torch.clamp(eigenvalues, min=0.0)
    L = V * torch.sqrt(eigenvalues)                       # L @ L^T = crlb_xyz
    z = torch.randn(3 * M, device=device, dtype=dtype)
    W = L @ z                                             # (3M,) ~ N(0, crlb_xyz)

    # Step 5: UGIA-F sample.
    xyzF = xyz + W.view(M, 3)
    return xyzF, F_xyz, crlb_xyz


# ============================================================
# Batched EM localization (3D astigmatic)
# ============================================================

def astig3d_em_batch_torch(
    frames:    torch.Tensor,
    xyz_init:  torch.Tensor,
    system:    SystemConfig,
    GN:        str = 'N',
    n_iter:    int = 10,
    dtype:     torch.dtype = torch.float32,
    z_clamp:   Optional[tuple] = None,
) -> torch.Tensor:
    """Batched 3D astigmatic EM localization with a fixed iteration count.

    The 3D analog of ``gauss2d_em_batch_torch``. Processes all B frames at
    once with a leading batch axis and a fixed ``n_iter`` outer EM loop, so
    there is no data-dependent stopping and the whole batch runs on the GPU
    in parallel. The only structural change from the 2D batch EM is the third
    coordinate z, whose gradient and single-emitter Fisher information come
    from the depth dependence of the astigmatic widths sigma_x(z), sigma_y(z),
    through d(sigma^2)/dz and the same closed-form (x^2-1)e^{-x^2/2} pixel
    integral used in ``astig3d_fisher_torch``. The M-step is a 3x3 noise-free
    single-emitter Newton solve, GN='N', or a diagonal-Fisher preconditioned
    gradient ascent, otherwise.

    Parameters
    ----------
    frames : (B, Ky, Kx) tensor
        Raw observed frames. The Gaussian-to-Poisson correction is applied
        inside this function.
    xyz_init : (B, M, 3) tensor
        Initial emitter positions (x, y, z) in nm for each of the B frames.
    system : SystemConfig
        Provides the Astigmatic3DPSF, Camera, Noise, and emitter intensities.
        ``system.frames`` is ignored; pass ``frames`` directly.
    GN : str
        'N' -> 3x3 Newton single-emitter M-step (default); else gradient.
    n_iter : int
        Number of EM outer iterations.
    dtype : torch.dtype
        float32 default; pass float64 for high photon counts.
    z_clamp : (zmin, zmax) or None
        If given, z is clamped to [zmin, zmax] each iteration, keeping
        emitters inside the calibrated depth range. x and y are always
        clamped to the frame.

    Returns
    -------
    xyz : (B, M, 3) tensor of estimated emitter positions in nm.
    """
    psf    = system.psf
    camera = system.camera
    noise  = system.noise
    if psf.psf_type != "astigmatic3d":
        raise ValueError("astig3d_em_batch_torch requires Astigmatic3DPSF.")

    c, d     = psf.c, psf.d
    sx0, sy0 = psf.sigmax0, psf.sigmay0
    Ax, Bx   = psf.Ax, psf.Bx
    Ay, By   = psf.Ay, psf.By

    Kx, Ky   = camera.Kx, camera.Ky
    Dx, Dy   = camera.Dx, camera.Dy
    Dt       = camera.Dt
    DxDyDt   = Dx * Dy * Dt

    device   = frames.device
    B        = frames.shape[0]
    M        = xyz_init.shape[1]

    b  = noise.b.to(device=device, dtype=dtype)    # (Ky, Kx)
    G  = noise.G.to(device=device, dtype=dtype)
    mu = noise.mu.to(device=device, dtype=dtype)

    # Gaussian-to-Poisson correction (all B frames at once).
    U = frames.to(dtype=dtype) + DxDyDt * (G - mu)             # (B, Ky, Kx)

    Im        = system.emitter.Im[0].to(device=device, dtype=dtype)   # (M,)
    Im_v      = Im.view(1, M, 1)                                # (1, M, 1)
    DtDxDyI_v = (DxDyDt * Im).view(1, M)                        # (1, M)

    alpha_newton   = 1.0
    alpha_gradient = 1.0;  iter_no_grad = 3
    sqrt2pi = math.sqrt(2.0 * math.pi)

    kx_v = torch.arange(Kx, device=device, dtype=dtype).view(1, 1, Kx)
    ky_v = torch.arange(Ky, device=device, dtype=dtype).view(1, 1, Ky)

    xyz0 = xyz_init.to(device=device, dtype=dtype).clone()     # (B, M, 3)
    xyz  = xyz0.clone()

    def _psf_and_grads(x0v, y0v, z0):
        # x0v, y0v : (B, M, 1);  z0 : (B, M)
        zx = (z0 + c) / d
        zy = (z0 - c) / d
        sx = sx0 * torch.sqrt((1.0 + zx**2 + Ax*zx**3 + Bx*zx**4).clamp_min(1e-12))
        sy = sy0 * torch.sqrt((1.0 + zy**2 + Ay*zy**3 + By*zy**4).clamp_min(1e-12))
        dsx2 = sx0**2 * (2*zx + 3*Ax*zx**2 + 4*Bx*zx**3) / d
        dsy2 = sy0**2 * (2*zy + 3*Ay*zy**2 + 4*By*zy**3) / d
        sxu = sx.unsqueeze(2); syu = sy.unsqueeze(2)            # (B, M, 1)

        Dxk1 = (Dx*(kx_v + 1.0) - x0v) / sxu                    # (B, M, Kx)
        Dxk0 = (Dx* kx_v        - x0v) / sxu
        Dyk1 = (Dy*(ky_v + 1.0) - y0v) / syu                    # (B, M, Ky)
        Dyk0 = (Dy* ky_v        - y0v) / syu

        qx = (Qfunc_torch(-Dxk1) - Qfunc_torch(-Dxk0)) / Dx
        qy = (Qfunc_torch(-Dyk1) - Qfunc_torch(-Dyk0)) / Dy
        DqxDx = (torch.exp(-0.5*Dxk0**2) - torch.exp(-0.5*Dxk1**2)) / (Dx*sqrt2pi*sxu)
        DqyDy = (torch.exp(-0.5*Dyk0**2) - torch.exp(-0.5*Dyk1**2)) / (Dy*sqrt2pi*syu)
        iF0x  = Dxk0*torch.exp(-0.5*Dxk0**2) - Dxk1*torch.exp(-0.5*Dxk1**2)
        iF0y  = Dyk0*torch.exp(-0.5*Dyk0**2) - Dyk1*torch.exp(-0.5*Dyk1**2)
        DqxDz = (dsx2.unsqueeze(2) / (2*Dx*sxu**2*sqrt2pi)) * iF0x
        DqyDz = (dsy2.unsqueeze(2) / (2*Dy*syu**2*sqrt2pi)) * iF0y

        q    = qy.unsqueeze(3)    * qx.unsqueeze(2)             # (B, M, Ky, Kx)
        dqdx = qy.unsqueeze(3)    * DqxDx.unsqueeze(2)
        dqdy = DqyDy.unsqueeze(3) * qx.unsqueeze(2)
        dqdz = DqxDz.unsqueeze(2) * qy.unsqueeze(3) + qx.unsqueeze(2) * DqyDz.unsqueeze(3)
        return q, dqdx, dqdy, dqdz

    for _ in range(n_iter):

        x0v = xyz0[:, :, 0].unsqueeze(2)                       # (B, M, 1)
        y0v = xyz0[:, :, 1].unsqueeze(2)
        z0  = xyz0[:, :, 2]                                    # (B, M)
        q, dqdx, dqdy, dqdz = _psf_and_grads(x0v, y0v, z0)

        sm = DxDyDt * Im_v.unsqueeze(3) * q                    # (B, M, Ky, Kx)
        u  = torch.clamp(sm.sum(dim=1) + DxDyDt * (b + G).unsqueeze(0), min=1e-10)
        Sm = sm * U.unsqueeze(1) / u.unsqueeze(1)              # E-step (B, M, Ky, Kx)

        s_all = torch.clamp(DtDxDyI_v.view(1, M, 1, 1) * q, min=1e-10)
        Ss1   = Sm / s_all - 1.0
        Lx = DtDxDyI_v * (Ss1 * dqdx).sum(dim=(2, 3))          # (B, M)
        Ly = DtDxDyI_v * (Ss1 * dqdy).sum(dim=(2, 3))
        Lz = DtDxDyI_v * (Ss1 * dqdz).sum(dim=(2, 3))
        Qrp = torch.clamp(q, min=1e-10)

        if GN == 'N':
            Fxx = DtDxDyI_v * (dqdx*dqdx / Qrp).sum(dim=(2, 3))   # (B, M)
            Fxy = DtDxDyI_v * (dqdx*dqdy / Qrp).sum(dim=(2, 3))
            Fxz = DtDxDyI_v * (dqdx*dqdz / Qrp).sum(dim=(2, 3))
            Fyy = DtDxDyI_v * (dqdy*dqdy / Qrp).sum(dim=(2, 3))
            Fyz = DtDxDyI_v * (dqdy*dqdz / Qrp).sum(dim=(2, 3))
            Fzz = DtDxDyI_v * (dqdz*dqdz / Qrp).sum(dim=(2, 3))

            Fmat = torch.stack([
                torch.stack([Fxx, Fxy, Fxz], dim=-1),
                torch.stack([Fxy, Fyy, Fyz], dim=-1),
                torch.stack([Fxz, Fyz, Fzz], dim=-1),
            ], dim=-2)                                          # (B, M, 3, 3)
            Lvec = torch.stack([Lx, Ly, Lz], dim=-1).unsqueeze(-1)   # (B, M, 3, 1)

            # Ridge for invertibility at the q->0 / weak-axial region.
            eye = torch.eye(3, device=device, dtype=dtype)
            ridge = 1e-6 * Fmat.diagonal(dim1=-2, dim2=-1).abs().amax(dim=-1, keepdim=True).unsqueeze(-1)
            Fmat = Fmat + (ridge + 1e-30) * eye

            delta = torch.linalg.solve(Fmat, Lvec).squeeze(-1)      # (B, M, 3)
            # Trust region: lateral step <= max width, axial step <= depth d.
            tr_xy = max(sx0, sy0)
            delta[:, :, 0] = (alpha_newton * delta[:, :, 0]).clamp(-tr_xy, tr_xy)
            delta[:, :, 1] = (alpha_newton * delta[:, :, 1]).clamp(-tr_xy, tr_xy)
            delta[:, :, 2] = (alpha_newton * delta[:, :, 2]).clamp(-d, d)
            xyz = xyz0 + delta
        else:
            xg = x0v.clone(); yg = y0v.clone(); zg = z0.clone()
            for _g in range(iter_no_grad):
                Fxxg = DtDxDyI_v * (dqdx*dqdx / Qrp).sum(dim=(2, 3))
                Fyyg = DtDxDyI_v * (dqdy*dqdy / Qrp).sum(dim=(2, 3))
                Fzzg = DtDxDyI_v * (dqdz*dqdz / Qrp).sum(dim=(2, 3))
                tr_xy = max(sx0, sy0)
                sxs = (alpha_gradient * Lx / (Fxxg + 1e-30)).clamp(-tr_xy, tr_xy)
                sys = (alpha_gradient * Ly / (Fyyg + 1e-30)).clamp(-tr_xy, tr_xy)
                szs = (alpha_gradient * Lz / (Fzzg + 1e-30)).clamp(-d, d)
                xg = (xg + sxs.unsqueeze(2)).clamp(0.0, Kx * Dx)
                yg = (yg + sys.unsqueeze(2)).clamp(0.0, Ky * Dy)
                zg = zg + szs
                if z_clamp is not None:
                    zg = zg.clamp(z_clamp[0], z_clamp[1])
                if _g < iter_no_grad - 1:
                    q, dqdx, dqdy, dqdz = _psf_and_grads(xg, yg, zg)
                    s_all = torch.clamp(DtDxDyI_v.view(1, M, 1, 1) * q, min=1e-10)
                    Ss1 = Sm / s_all - 1.0
                    Lx = DtDxDyI_v * (Ss1 * dqdx).sum(dim=(2, 3))
                    Ly = DtDxDyI_v * (Ss1 * dqdy).sum(dim=(2, 3))
                    Lz = DtDxDyI_v * (Ss1 * dqdz).sum(dim=(2, 3))
                    Qrp = torch.clamp(q, min=1e-10)
            xyz = torch.stack([xg.squeeze(2), yg.squeeze(2), zg], dim=-1)

        # Keep emitters on the frame and inside the calibrated depth.
        xyz[:, :, 0] = xyz[:, :, 0].clamp(0.0, Kx * Dx)
        xyz[:, :, 1] = xyz[:, :, 1].clamp(0.0, Ky * Dy)
        if z_clamp is not None:
            xyz[:, :, 2] = xyz[:, :, 2].clamp(z_clamp[0], z_clamp[1])

        xyz0 = xyz.clone()

    return xyz   # (B, M, 3)


# ============================================================
# Information-sufficient SNR (3D astigmatic, per emitter, per frame)
# ============================================================

def astig3d_snr_torch(
    system: SystemConfig,
    chunk_size: Optional[int] = None,
) -> torch.Tensor:
    """Per-emitter information-sufficient SNR for the 3D astigmatic PSF.

    The 3D analog of ``gauss2d_snr_torch``. Returns ``snr_db`` of shape
    (N, M), the 80%-power-effective SNR of Sun (2020) extended by M. Sun and
    Y. Sun (2022) to the astigmatic PSF, with unequal emitter intensities and
    per-pixel non-uniform noise:

        n_bar_m = sum_k q_m(k) (b_k + G_k) / sum_k q_m(k)            (Eq. 27)
        gamma_m = I_m / n_bar_m
        eta_m   = eta_s / (sigma_x(z_m) sigma_y(z_m))                (Eq. 30)
        SNR_m   = 10 log10(gamma_m) - 10 log10(sigma_x sigma_y) - 11.02  (Eq. 31)

    Here q_m(k) is the pixel-integrated astigmatic PSF, with the depth-
    dependent widths sigma_x(z_m), sigma_y(z_m) of Appendix C entering the x
    and y integrals separately; b_k, G_k are the background and readout
    variance densities in photons/(s*nm^2); I_m is the intensity in
    photons/s. The figure-of-merit constant is eta_s = -rho/[2 pi ln(1-rho)]
    = 0.0791 at rho = 0.80, so 10 log10(eta_s) = -11.02 dB. Inactive emitters
    (Im == 0) return -inf.

    Parameters
    ----------
    system : SystemConfig
        Must contain an Astigmatic3DPSF, a 3-D EmitterData, a Camera, and a
        Noise.
    chunk_size : int, optional
        Process the N frames chunk_size at a time. None does all N at once.

    Returns
    -------
    snr_db : (N, M) tensor — per-emitter SNR in dB. Per-frame range and mean
        over active emitters via ``snr_db.amin(1)``, ``amax(1)``, ``mean(1)``.

    References
    ----------
    Sun, Y. Opt. Lett. 45(21), 6102-6105 (2020).
    M. Sun and Y. Sun, Electron. Lett. 58(2), 58-60 (2022), Eq. (9).
    """
    emitter = system.emitter
    psf = system.psf
    camera = system.camera
    noise = system.noise
    if psf.psf_type != "astigmatic3d":
        raise ValueError(
            f"astig3d_snr_torch requires Astigmatic3DPSF; got '{psf.psf_type}'."
        )

    xyz = emitter.positions          # (N, M, 3)
    Im = emitter.Im                  # (N, M)
    N, M, _ = xyz.shape
    device, dtype = xyz.device, xyz.dtype

    Kx, Ky = camera.Kx, camera.Ky
    Dx, Dy = camera.Dx, camera.Dy

    Im = Im.to(device=device, dtype=dtype)
    n_density = (noise.b + noise.G).to(device=device, dtype=dtype)   # (Ky, Kx)
    sqrt2 = math.sqrt(2.0)

    kx = torch.arange(Kx, device=device, dtype=dtype).view(1, 1, 1, Kx)
    ky = torch.arange(Ky, device=device, dtype=dtype).view(1, 1, Ky, 1)
    Dxk0, Dxk1 = Dx * kx, Dx * (kx + 1)
    Dyk0, Dyk1 = Dy * ky, Dy * (ky + 1)

    rho = 0.80
    eta_s = -rho / (2.0 * math.pi * math.log(1.0 - rho))   # 0.07911...

    snr_db = torch.empty((N, M), device=device, dtype=dtype)
    chunk = N if chunk_size is None else int(chunk_size)
    if chunk <= 0:
        raise ValueError("chunk_size must be a positive integer or None.")

    for n0 in range(0, N, chunk):
        n1 = min(n0 + chunk, N)

        x0 = xyz[n0:n1, :, 0].view(-1, M, 1, 1)
        y0 = xyz[n0:n1, :, 1].view(-1, M, 1, 1)
        z  = xyz[n0:n1, :, 2]                              # (Nc, M)
        sx, sy = astig_sigma_xy(z, psf)                    # (Nc, M)
        sxu = sx.view(-1, M, 1, 1); syu = sy.view(-1, M, 1, 1)

        qx = 0.5 * (torch.erf((Dxk1 - x0) / (sxu * sqrt2))
                    - torch.erf((Dxk0 - x0) / (sxu * sqrt2)))   # (Nc, M, 1, Kx)
        qy = 0.5 * (torch.erf((Dyk1 - y0) / (syu * sqrt2))
                    - torch.erf((Dyk0 - y0) / (syu * sqrt2)))   # (Nc, M, Ky, 1)
        q = qy * qx                                            # (Nc, M, Ky, Kx)

        num = (q * n_density).sum(dim=(-2, -1))
        den = q.sum(dim=(-2, -1))
        n_bar = num / den                                     # (Nc, M)

        gamma = Im[n0:n1] / n_bar                             # (Nc, M)
        eta = eta_s / (sx * sy)                               # (Nc, M)
        snr_db[n0:n1] = 10.0 * torch.log10(eta * gamma)

    active = Im > 0
    return torch.where(active, snr_db, torch.full_like(snr_db, float("-inf")))
