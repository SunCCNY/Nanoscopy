"""
config.py

Dataclass configuration for physical SMLM system parameters.

Conventions
-----------
- Arrays/tensors are torch.Tensor (PyTorch chosen for downstream NN/ML work
  and free GPU/autodiff support).
- Spatial units: nanometers (nm).
- Temporal units: seconds (s).
- Photon densities: photons / (s * nm^2).
- Image arrays use (Ky, Kx) shape — i.e., (rows, cols) = (height, width).
- Emitter data uses (N, M, ...) — N frames * M emitters.

Choice of device (CPU/GPU) and dtype (float32/float64) is left entirely
to the caller. Construct each tensor field on the desired device when
building the config; functions that consume the config are responsible
for any device alignment they need.

Reference
---------
Data-frame model: Sun (2013, JBO) and Sun & Guan (2021, JOSA A), with two
extensions adopted in this library:
  (1) the Poisson-background density `b` and the Gaussian-readout
      density/variance `mu`, `G` are spatially variant (per-pixel) but
      temporally invariant;
  (2) emitter locations may differ across frames (the within-frame model
      is unchanged: emitters are stationary during one exposure).
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
import math

import torch


# ============================================================
# Internal helpers
# ============================================================

def _promote_to_map(x, ref: torch.Tensor) -> torch.Tensor:
    """Promote a scalar or 0-D tensor to a 2-D map matching ``ref`` in
    shape, dtype, and device. A 2-D tensor is returned unchanged.

    Used so that uniform Gaussian readout (a common simplification or the
    EMCCD case) can be specified as a single number while spatially-varying
    readout (sCMOS) is specified as a (Ky, Kx) tensor.
    """
    if isinstance(x, torch.Tensor):
        if x.ndim == 0:
            return torch.full_like(ref, x.item())
        if x.ndim == 2:
            return x
        raise ValueError(
            f"Expected scalar, 0-D tensor, or 2-D tensor; got ndim={x.ndim}."
        )
    return torch.full_like(ref, float(x))


# ============================================================
# Emitter Data
# ============================================================

@dataclass
class EmitterData:
    """Emitter trajectories and intensities for an SMLM movie.

    Parameters
    ----------
    xy : (N, M, 2) tensor, optional
        2D emitter positions in nm.
    xyz : (N, M, 3) tensor, optional
        3D emitter positions in nm. Exactly one of `xy` or `xyz` must
        be provided.
    Im : (N, M) tensor
        Per-frame, per-emitter intensities in photons/s. ``Im[n, m] == 0``
        marks emitter ``m`` as inactive in frame ``n``.

    Notes
    -----
    A single frame is represented by N == 1. Within each frame, emitter
    locations are stationary (intra-frame motion blur is not modelled in
    this layer); across frames they may differ to support photoactivation
    kinetics or true emitter motion handled by upper layers.
    """

    xy: Optional[torch.Tensor] = None
    xyz: Optional[torch.Tensor] = None
    Im: Optional[torch.Tensor] = None

    def __post_init__(self):
        # Exactly one of xy or xyz
        if (self.xy is None) == (self.xyz is None):
            raise ValueError("Provide exactly one of xy or xyz, not both.")

        # Validate Im
        if not isinstance(self.Im, torch.Tensor) or self.Im.ndim != 2:
            raise ValueError("Im must be a 2-D tensor of shape (N, M).")

        # Validate positions
        if self.xy is not None:
            if (not isinstance(self.xy, torch.Tensor) or self.xy.ndim != 3
                    or self.xy.shape[-1] != 2):
                raise ValueError("xy must be a tensor of shape (N, M, 2).")
            N1, M1, _ = self.xy.shape
        else:
            if (not isinstance(self.xyz, torch.Tensor) or self.xyz.ndim != 3
                    or self.xyz.shape[-1] != 3):
                raise ValueError("xyz must be a tensor of shape (N, M, 3).")
            N1, M1, _ = self.xyz.shape

        # Cross-validate shapes
        N2, M2 = self.Im.shape
        if (N1, M1) != (N2, M2):
            raise ValueError(
                f"Position and Im must share (N, M); got positions={(N1, M1)}, "
                f"Im={(N2, M2)}."
            )

        # Physical sanity
        if torch.any(self.Im < 0):
            raise ValueError("Im must be non-negative.")

    # --- shape properties ---

    @property
    def N(self) -> int:
        return self.Im.shape[0]

    @property
    def M(self) -> int:
        return self.Im.shape[1]

    @property
    def dim(self) -> int:
        return 2 if self.xy is not None else 3

    @property
    def positions(self) -> torch.Tensor:
        """(N, M, 2) for 2-D, (N, M, 3) for 3-D."""
        return self.xy if self.xy is not None else self.xyz

    @property
    def active(self) -> torch.Tensor:
        """Boolean (N, M) mask: True where the emitter is on (Im > 0)."""
        return self.Im > 0


# ============================================================
# PSF Models
# ============================================================

@dataclass
class PSF:
    """Base class for PSF descriptors. Subclasses set ``psf_type`` and
    expose their parameters."""
    psf_type: str

    @property
    def dim(self) -> int:
        """Spatial dimensionality of the PSF (2 or 3)."""
        raise NotImplementedError


@dataclass
class Gaussian2DPSF(PSF):
    """2-D isotropic Gaussian PSF.

    Parameters
    ----------
    sigma : float
        Standard deviation in nm (same in x and y).
    """
    sigma: float

    def __init__(self, sigma: float):
        super().__init__(psf_type="gaussian2d")
        self.sigma = float(sigma)

    @property
    def dim(self) -> int:
        return 2


@dataclass
class Airy2DPSF(PSF):
    """2-D Airy PSF at the focal plane (circular aperture).

        q(v) = J1^2(alpha * |v - v_m|) / (pi * |v - v_m|^2),
        alpha = 2 * pi * na / lambda

    Parameters
    ----------
    na : float
        Numerical aperture (dimensionless).
    lambda_ : float
        Fluorescence emission wavelength in nm.

    See Also
    --------
    equivalent_gaussian_sigma : Gaussian SD that matches this Airy PSF
        in the sense of Sun (2013, JBO) Appendix B.
    """
    na: float
    lambda_: float

    def __init__(self, na: float, lambda_: float):
        super().__init__(psf_type="airy2d")
        self.na = float(na)
        self.lambda_ = float(lambda_)

    @property
    def dim(self) -> int:
        return 2

    @property
    def alpha(self) -> float:
        """alpha = 2 * pi * na / lambda  (units: 1/nm)."""
        return 2.0 * math.pi * self.na / self.lambda_

    @property
    def equivalent_gaussian_sigma(self) -> float:
        """Equivalent 2-D Gaussian SD in nm, ~= 1.3238 / alpha.

        From Sun (2013, JBO) Appendix B, integrating J1^2 up to its
        second zero (alpha * x ~= 7.016).
        """
        return 1.3238 / self.alpha


@dataclass
class Astigmatic3DPSF(PSF):
    """3-D astigmatic Gaussian PSF for cylindrical-lens 3-D nanoscopy.

        sigma_x(z) = sigmax0 * sqrt(1 + (z+c)^2/d^2 + Ax(z+c)^3/d^3 + Bx(z+c)^4/d^4)
        sigma_y(z) = sigmay0 * sqrt(1 + (z-c)^2/d^2 + Ay(z-c)^3/d^3 + By(z-c)^4/d^4)
    """
    c: float
    d: float
    sigmax0: float
    Ax: float
    Bx: float
    sigmay0: float
    Ay: float
    By: float

    def __init__(self, c, d, sigmax0, Ax, Bx, sigmay0, Ay, By):
        super().__init__(psf_type="astigmatic3d")
        self.c = float(c)
        self.d = float(d)
        self.sigmax0 = float(sigmax0)
        self.Ax = float(Ax)
        self.Bx = float(Bx)
        self.sigmay0 = float(sigmay0)
        self.Ay = float(Ay)
        self.By = float(By)

    @property
    def dim(self) -> int:
        return 3


# ============================================================
# Camera Geometry
# ============================================================

@dataclass
class Camera:
    """Camera geometry. Noise is handled by the Noise class.

    Parameters
    ----------
    Kx, Ky : int
        Frame width and height in pixels.
    Dx, Dy : float
        Pixel sizes in nm.
    Dt : float
        Frame exposure time in seconds.
    """

    Kx: int
    Ky: int
    Dx: float
    Dy: float
    Dt: float

    def __post_init__(self):
        if self.Kx <= 0 or self.Ky <= 0:
            raise ValueError("Kx and Ky must be positive integers.")
        if self.Dx <= 0 or self.Dy <= 0:
            raise ValueError("Dx and Dy must be positive.")
        if self.Dt <= 0:
            raise ValueError("Dt must be positive.")
        # normalize types
        self.Kx = int(self.Kx)
        self.Ky = int(self.Ky)
        self.Dx = float(self.Dx)
        self.Dy = float(self.Dy)
        self.Dt = float(self.Dt)

    @property
    def pixel_area(self) -> float:
        """Dx * Dy (nm^2)."""
        return self.Dx * self.Dy

    @property
    def Lx(self) -> float:
        """Field of view width in nm."""
        return self.Kx * self.Dx

    @property
    def Ly(self) -> float:
        """Field of view height in nm."""
        return self.Ky * self.Dy


# ============================================================
# Noise Model (Poisson + Gaussian, all per-pixel)
# ============================================================

@dataclass
class Noise:
    """Per-pixel noise model.

    Parameters
    ----------
    b : (Ky, Kx) tensor
        Spatially varying mean density of background Poisson noise
        (autofluorescence). Units: photons / (s * nm^2).
    mu : (Ky, Kx) tensor or scalar
        Per-pixel mean density of Gaussian readout noise. For sCMOS
        sensors this typically varies pixel-to-pixel; for EMCCDs (or as
        a simplification) a scalar may be passed and will be broadcast
        to a uniform map matching `b`. Units: photons / (s * nm^2).
    G : (Ky, Kx) tensor or scalar
        Per-pixel variance density of Gaussian readout noise. Same
        per-pixel-vs-uniform handling as `mu`. Units: photons / (s * nm^2).

    Notes
    -----
    All three maps share (Ky, Kx). The mean `mu` of Gaussian noise does
    not affect the Cramer-Rao bound (it can be subtracted from the data
    before estimation), but it is retained here so that simulated frames
    reproduce raw pixel values realistically.
    """

    b: torch.Tensor
    mu: torch.Tensor
    G: torch.Tensor

    def __post_init__(self):
        if not isinstance(self.b, torch.Tensor) or self.b.ndim != 2:
            raise ValueError("b must be a 2-D tensor of shape (Ky, Kx).")

        # Promote scalar mu/G (uniform readout) to a 2-D map matching b.
        self.mu = _promote_to_map(self.mu, self.b)
        self.G = _promote_to_map(self.G, self.b)

        if self.mu.shape != self.b.shape:
            raise ValueError(
                f"mu shape {tuple(self.mu.shape)} != b shape {tuple(self.b.shape)}."
            )
        if self.G.shape != self.b.shape:
            raise ValueError(
                f"G shape {tuple(self.G.shape)} != b shape {tuple(self.b.shape)}."
            )

        if torch.any(self.b < 0):
            raise ValueError("b must be non-negative.")
        if torch.any(self.G < 0):
            raise ValueError("G (variance) must be non-negative.")

    # --- shape properties ---

    @property
    def Ky(self) -> int:
        return self.b.shape[0]

    @property
    def Kx(self) -> int:
        return self.b.shape[1]

    # --- per-pixel integrated quantities ---
    # Note: these return (Ky, Kx) tensors on the same device/dtype as the
    # noise tensors. Callers needing a different device should add their
    # own .to(device) — see SMLM_Lib functions for examples.

    def poisson_mean_per_pixel(self, camera: Camera) -> torch.Tensor:
        """Dt * Dx * Dy * b -- mean Poisson background photons per pixel
        (Ky, Kx)."""
        return self.b * (camera.Dt * camera.pixel_area)

    def gaussian_mean_per_pixel(self, camera: Camera) -> torch.Tensor:
        """Dt * Dx * Dy * mu -- mean of Gaussian readout per pixel (Ky, Kx)."""
        return self.mu * (camera.Dt * camera.pixel_area)

    def gaussian_var_per_pixel(self, camera: Camera) -> torch.Tensor:
        """Dt * Dx * Dy * G -- variance of Gaussian readout per pixel
        (Ky, Kx)."""
        return self.G * (camera.Dt * camera.pixel_area)


# ============================================================
# System Configuration
# ============================================================

@dataclass
class SystemConfig:
    """Bundles emitters, PSF, camera, and noise into a single config.

    Parameters
    ----------
    emitter, psf, camera, noise : the four required components.
    frames : (N, Ky, Kx) tensor, optional
        Acquired or simulated movie data. A single frame is represented
        by N == 1.
    """

    emitter: EmitterData
    psf: PSF
    camera: Camera
    noise: Noise

    frames: Optional[torch.Tensor] = None

    def __post_init__(self):
        if self.emitter is None:
            raise ValueError("SystemConfig requires an emitter.")
        if self.psf is None:
            raise ValueError("SystemConfig requires a PSF model.")
        if self.camera is None:
            raise ValueError("SystemConfig requires a camera model.")
        if self.noise is None:
            raise ValueError("SystemConfig requires a Noise model.")

        if not hasattr(self.psf, "psf_type"):
            raise ValueError("PSF must define psf_type.")

        # Cross-validate: noise maps must match camera (Ky, Kx).
        if (self.noise.Ky, self.noise.Kx) != (self.camera.Ky, self.camera.Kx):
            raise ValueError(
                f"Noise shape (Ky, Kx) = ({self.noise.Ky}, {self.noise.Kx}) "
                f"does not match camera ({self.camera.Ky}, {self.camera.Kx})."
            )

        # PSF dimensionality must match emitter dimensionality.
        if self.psf.dim != self.emitter.dim:
            raise ValueError(
                f"PSF '{self.psf.psf_type}' is {self.psf.dim}-D but emitter "
                f"is {self.emitter.dim}-D."
            )

        # Validate frames shape if provided.
        if self.frames is not None:
            if not isinstance(self.frames, torch.Tensor) or self.frames.ndim != 3:
                raise ValueError("frames must be a 3-D tensor (N, Ky, Kx).")
            N, Ky, Kx = self.frames.shape
            if (Ky, Kx) != (self.camera.Ky, self.camera.Kx):
                raise ValueError(
                    f"frames image size ({Ky}, {Kx}) does not match camera "
                    f"({self.camera.Ky}, {self.camera.Kx})."
                )
            if N != self.emitter.N:
                raise ValueError(
                    f"frames N={N} does not match emitter N={self.emitter.N}."
                )
