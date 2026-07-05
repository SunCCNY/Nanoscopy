"""
SMLM_Lib: A Python Library for Single-Molecule Localization Microscopy (SMLM)
"""

# ------------------------------------------------------------
# Configuration classes
# ------------------------------------------------------------
from .config import (
    EmitterData,
    PSF,                  # base class — useful for subclassing new PSF types
    Gaussian2DPSF,
    Airy2DPSF,            # was missing
    Astigmatic3DPSF,
    Camera,
    Noise,
    SystemConfig,
)

# ------------------------------------------------------------
# Emitter activation models
# ------------------------------------------------------------
from .activation import (
    emActMarkovContinue,
    emActMarkovCycle,
    activation_statistics_continue,
    activation_statistics_cycle,
)

# ------------------------------------------------------------
# Astig3D forward model + estimators
# ------------------------------------------------------------
from .Astig3D import (
    astig3d_frame_torch,
    astig_sigma_xy,
    astig3d_fisher_torch,
    astig3d_ugia_f_torch,
    astig3d_em_batch_torch,
    astig3d_snr_torch,
)

# ------------------------------------------------------------
# Gauss2D forward model + estimators
# ------------------------------------------------------------
from .Gauss2D import (
    Qfunc_torch,
    gauss2d_frame_torch,
    gauss2d_fisher_torch,
    gauss2d_ugia_f_torch,
    gauss2d_ugia_f_batch_torch,
    gauss2d_seml_newton,
    gauss2d_seml_gradient,
    gauss2d_em_torch,
    gauss2d_em_batch_torch,
    gauss2d_snr_torch,
)

# ------------------------------------------------------------
# Metrics
# ------------------------------------------------------------
from .metrics import (
    rmsmd,
    rmsmd_per_frame,
    partition_x,
    rmse_p,
    rmsmd_p,
)

# ------------------------------------------------------------
# Emitter position utilities
# ------------------------------------------------------------
from .position import (
    sample_emitters_batch,
    circle_emitters,
    background_cloud,
    sample_delta_separated,
    sample_delta_separated_3d,
)

# ------------------------------------------------------------
# Shared training utilities (batch generator + Hungarian loss)
# ------------------------------------------------------------
from .train_utils import (
    generate_batch,
    hungarian_loss,
)

# ------------------------------------------------------------
# APL models (formerly "KAN"): linear mixing + learnable
# piecewise-linear activation  (Agostinelli/Scardapane lineage)
# ------------------------------------------------------------
from .APL import (
    APLLayer,
    APLLocalizer,
    CNNAPLLocalizer,
)

# ------------------------------------------------------------
# True KAN models (Liu et al. 2024): B-spline edge functions
# ------------------------------------------------------------
from .KAN import (
    KANLayer,
    KANLocalizer,
    CNNKANLocalizer,
)

# ------------------------------------------------------------
# Visualization utilities
# ------------------------------------------------------------
from .visualize import (
    show_frame_with_localizations,
    print_device_info,
)

# ------------------------------------------------------------
# ViT localizer (Hungarian-based)
# ------------------------------------------------------------
from .ViT import (
    ViTLocalizer,
)

# ------------------------------------------------------------
# ViT-PF localizer (Permutation-Free, subpixel classification)
# ------------------------------------------------------------
from .ViT_PF import (
    ViTPFLocalizer,
    generate_batch_pf,
)
