"""
Emitter position generation utilities for SMLM_Lib.
All coordinates are in nanometers unless otherwise stated.
"""

import math

import torch
import numpy as np

from .config import Camera


# ============================================================
# 1. Random emitter sampling inside the frame (nm)
# ============================================================

def sample_emitters_batch(batch_size, M, camera: Camera, margin_px=2):
    """
    Sample random emitter positions inside the frame, in nanometers.

    Args:
        batch_size : number of frames
        M          : number of emitters per frame
        camera     : Camera object (provides Kx, Ky, Dx, Dy)
        margin_px  : avoid edges by this many pixels

    Returns:
        xy : (batch_size, M, 2) tensor in nanometers
    """

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    Kx, Ky = camera.Kx, camera.Ky
    Dx, Dy = camera.Dx, camera.Dy

    # Allowed pixel indices
    allowed_x = torch.arange(margin_px, Kx - margin_px, device=device)
    allowed_y = torch.arange(margin_px, Ky - margin_px, device=device)

    # Random integer pixel positions
    px = allowed_x[torch.randint(0, len(allowed_x), (batch_size, M), device=device)]
    py = allowed_y[torch.randint(0, len(allowed_y), (batch_size, M), device=device)]

    # Random subpixel offsets
    rx = torch.rand(batch_size, M, device=device)
    ry = torch.rand(batch_size, M, device=device)

    # Convert to nanometers
    x_nm = (px + rx) * Dx
    y_nm = (py + ry) * Dy

    return torch.stack([x_nm, y_nm], dim=-1)   # (batch_size, M, 2)


# ============================================================
# 2. Circle emitter generation (nm)
# ============================================================

def circle_emitters(M, d_nm, camera: Camera, xy0_nm=(2, 5), th0=np.pi/3):
    """
    Generate M emitter positions on a circle of radius d_nm (nm).

    The circle is centered at xy0_nm (nm), then shifted to the
    center of the Kx × Ky frame.

    Returns:
        xy_nm : (M, 2) numpy array in nanometers
    """

    theta = np.arange(M) * 2 * np.pi / M + th0
    xy_nm = np.column_stack([
        d_nm * np.cos(theta) + xy0_nm[0],
        d_nm * np.sin(theta) + xy0_nm[1]
    ]).astype(np.float32)

    # Shift to center of frame
    shift_nm = np.array([
        (camera.Kx / 2) * camera.Dx,
        (camera.Ky / 2) * camera.Dy
    ], dtype=np.float32)

    return xy_nm + shift_nm

# ============================================================
# 3. Smooth background density field b_k (autofluorescence cloud)
# ============================================================

def background_cloud(camera: Camera, b_min, b_max, corr_px=10.0,
                     seed=None, device=None, dtype=torch.float32):
    """Smooth per-pixel background density field b_k for a frame.

    Models the slowly varying autofluorescence and out-of-focus haze as a
    stationary isotropic Gaussian random field: white noise on the
    (Ky, Kx) grid is low-pass filtered by a Gaussian of standard deviation
    ``corr_px`` pixels (the spatial correlation length, the cloud scale),
    then affinely mapped so the field spans [b_min, b_max] with mean near
    their midpoint. Larger ``corr_px`` gives a smoother, larger-scale cloud.

    Unlike the per-pixel-independent Gaussian readout G_k (sCMOS), the
    Poisson background varies slowly in space, so this field is smooth. It
    is one static map shared by all frames of a movie; the per-frame Poisson
    counts are sampled from it downstream.

    Args:
        camera   : Camera (provides Kx, Ky).
        b_min, b_max : range of the background density, photons/(s*nm^2).
        corr_px  : spatial correlation length in pixels.
        seed     : int or None. If given, the field is reproducible and
                   independent of the current global RNG state.
        device, dtype : placement of the returned tensor.

    Returns:
        b : (Ky, Kx) tensor — the background density map b_k.
    """
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    Kx, Ky = camera.Kx, camera.Ky

    gen = None
    if seed is not None:
        gen = torch.Generator(device="cpu").manual_seed(int(seed))

    # White noise, then Gaussian low-pass in the Fourier domain. Computed
    # in float64 on CPU for a device-independent, reproducible field.
    white = torch.randn(Ky, Kx, generator=gen, dtype=torch.float64)
    fy = torch.fft.fftfreq(Ky, dtype=torch.float64).view(Ky, 1)
    fx = torch.fft.fftfreq(Kx, dtype=torch.float64).view(1, Kx)
    f2 = fy * fy + fx * fx
    H = torch.exp(-2.0 * (math.pi ** 2) * (float(corr_px) ** 2) * f2)
    field = torch.fft.ifft2(torch.fft.fft2(white) * H).real

    field = (field - field.min()) / (field.max() - field.min())   # -> [0, 1]
    b = b_min + (b_max - b_min) * field
    return b.to(device=device, dtype=dtype)


# ============================================================
# 4. Delta-separated random emitter sampling (nm)
# ============================================================

def sample_delta_separated(M, camera: Camera, delta, margin_px=2, seed=0,
                           device=None, max_tries=2000000):
    """Rejection-sample M emitter positions (nm) with all pairwise distances
    at least ``delta`` nm.

    Candidates are drawn uniformly inside the frame (minus a ``margin_px``
    border) and accepted only if they are at least ``delta`` from every
    previously accepted emitter, so the returned configuration is
    delta-separated and resolvable while remaining dense. Reproducible from
    ``seed`` and independent of the global RNG state.

    Args:
        M         : number of emitters.
        camera    : Camera (provides Kx, Ky, Dx, Dy).
        delta     : minimum pairwise separation in nm.
        margin_px : keep emitters this many pixels from the frame edge.
        seed      : int, reproducible configuration.
        device    : optional device for the returned tensor.
        max_tries : safety cap on rejection attempts.

    Returns:
        xy : (M, 2) tensor in nanometers.

    Raises:
        RuntimeError if M emitters cannot be placed within max_tries (delta or
        M too large for the field of view).
    """
    g = torch.Generator().manual_seed(int(seed))
    x_lo, x_hi = margin_px * camera.Dx, (camera.Kx - margin_px) * camera.Dx
    y_lo, y_hi = margin_px * camera.Dy, (camera.Ky - margin_px) * camera.Dy
    pts = torch.empty(M, 2)
    n, tries, d2 = 0, 0, float(delta) * float(delta)
    while n < M and tries < max_tries:
        tries += 1
        c = torch.rand(2, generator=g)
        xy = torch.tensor([x_lo + c[0] * (x_hi - x_lo),
                           y_lo + c[1] * (y_hi - y_lo)])
        if n == 0 or torch.min(((pts[:n] - xy) ** 2).sum(1)).item() >= d2:
            pts[n] = xy
            n += 1
    if n < M:
        raise RuntimeError(
            f"Placed only {n}/{M} emitters with delta={delta} nm. "
            f"Lower delta or M, or enlarge the FOV.")
    return pts if device is None else pts.to(device)


# ============================================================
# 5. Delta-separated 3D emitter sampling (nm)
# ============================================================

def sample_delta_separated_3d(M, camera: Camera, delta, Lz, margin_px=2, seed=0,
                              device=None, max_tries=2000000,
                              delta_xy=None, delta_z=None):
    """Rejection-sample M 3D emitter positions (nm) with an anisotropic minimum
    separation.

    The 3D analog of ``sample_delta_separated``. Candidates are drawn uniformly
    inside the cuboid
        [margin_px*Dx, (Kx-margin_px)*Dx] x [margin_px*Dy, (Ky-margin_px)*Dy]
        x [-Lz, Lz].
    A candidate is accepted only if it lies outside the resolution ellipsoid of
    every previously accepted emitter,
        (dx^2 + dy^2) / delta_xy^2  +  dz^2 / delta_z^2  >=  1,
    so two emitters may be close laterally if well separated axially and vice
    versa. This matches the anisotropic resolution of the astigmatic PSF, whose
    axial information is much weaker than its lateral information, so a larger
    delta_z (and smaller delta_xy) keeps the configuration resolvable while
    remaining dense. With delta_xy = delta_z = delta the criterion reduces to the
    isotropic Euclidean test, all pairwise 3D distances at least ``delta``.
    Reproducible from ``seed`` and independent of the global RNG state.

    Args:
        M         : number of emitters.
        camera    : Camera (provides Kx, Ky, Dx, Dy).
        delta     : isotropic minimum separation in nm; used for whichever of
                    delta_xy, delta_z is left as None.
        Lz        : half axial range; z is drawn in [-Lz, Lz] nm.
        margin_px : keep emitters this many pixels from the lateral edge.
        seed      : int, reproducible configuration.
        device    : optional device for the returned tensor.
        max_tries : safety cap on rejection attempts.
        delta_xy  : lateral semi-axis of the resolution ellipsoid in nm
                    (defaults to ``delta``).
        delta_z   : axial semi-axis of the resolution ellipsoid in nm
                    (defaults to ``delta``).

    Returns:
        xyz : (M, 3) tensor in nanometers.

    Raises:
        RuntimeError if M emitters cannot be placed within max_tries (the
        separations or M are too large for the field of view).
    """
    delta_xy = float(delta if delta_xy is None else delta_xy)
    delta_z = float(delta if delta_z is None else delta_z)
    inv_xy2, inv_z2 = 1.0 / (delta_xy * delta_xy), 1.0 / (delta_z * delta_z)
    g = torch.Generator().manual_seed(int(seed))
    x_lo, x_hi = margin_px * camera.Dx, (camera.Kx - margin_px) * camera.Dx
    y_lo, y_hi = margin_px * camera.Dy, (camera.Ky - margin_px) * camera.Dy
    pts = torch.empty(M, 3)
    n, tries = 0, 0
    while n < M and tries < max_tries:
        tries += 1
        c = torch.rand(3, generator=g)
        xyz = torch.tensor([x_lo + c[0] * (x_hi - x_lo),
                            y_lo + c[1] * (y_hi - y_lo),
                            -Lz + c[2] * (2.0 * Lz)])
        if n == 0:
            pts[0] = xyz
            n = 1
            continue
        diff = pts[:n] - xyz
        sep = (diff[:, 0] ** 2 + diff[:, 1] ** 2) * inv_xy2 + diff[:, 2] ** 2 * inv_z2
        if torch.min(sep).item() >= 1.0:
            pts[n] = xyz
            n += 1
    if n < M:
        raise RuntimeError(
            f"Placed only {n}/{M} emitters with delta_xy={delta_xy} nm, "
            f"delta_z={delta_z} nm in 3D. Lower the separations or M, or "
            f"enlarge the FOV.")
    return pts if device is None else pts.to(device)
