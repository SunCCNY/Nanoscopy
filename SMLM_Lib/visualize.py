"""
visualize.py

Visualization utilities for SMLM localization results.
"""

import torch
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import zoom

from .config import Camera, Gaussian2DPSF, Noise, EmitterData, SystemConfig
from .Gauss2D import gauss2d_frame_torch
from .metrics import rmsmd_per_frame


def print_device_info():
    """Print GPU/CPU device status — call at start of every script."""
    if torch.cuda.is_available():
        gpu_name = torch.cuda.get_device_name(0)
        gpu_mem  = torch.cuda.get_device_properties(0).total_memory / 1024**3
        print(f"  Device : GPU — {gpu_name} ({gpu_mem:.1f} GB)")
    else:
        print("  Device : CPU  (no CUDA GPU available)")
    print(f"  PyTorch: {torch.__version__}")


def show_one_frame(frame, xy_nm, camera: Camera, imagej_coords=True):
    """Display a single frame with true emitter locations.

    Parameters
    ----------
    frame         : (Ky, Kx) numpy array or tensor
    xy_nm         : (M, 2) array/tensor — emitter positions in nm
    camera        : Camera
    imagej_coords : bool — if True use upper origin (ImageJ convention)
    """
    if isinstance(frame, torch.Tensor):
        frame = frame.cpu().numpy()
    if isinstance(xy_nm, torch.Tensor):
        xy_nm = xy_nm.cpu().numpy()

    xy_pix = xy_nm / np.array([camera.Dx, camera.Dy], dtype=np.float32) - 0.5
    origin = 'upper' if imagej_coords else 'lower'

    plt.figure(figsize=(5, 5))
    plt.imshow(frame, cmap='gray', origin=origin)
    plt.scatter(xy_pix[:, 0], xy_pix[:, 1], c='red', s=40)
    plt.title("Example Frame with True Emitter Locations")
    plt.tight_layout()
    plt.show()


def show_frame_with_localizations(
    model:         torch.nn.Module,
    M:             int,
    camera:        Camera,
    psf:           Gaussian2DPSF,
    noise:         Noise,
    true_xy_nm:    torch.Tensor,
    Im0:           float = 300000.0,
    imagej_coords: bool  = True,
):
    """Generate a frame from true positions, run the model, and plot results.

    Parameters
    ----------
    model       : trained localizer (APL / KAN / ViT) returning (B, M, 2)
    M           : number of emitters
    camera      : Camera
    psf         : Gaussian2DPSF
    noise       : Noise
    true_xy_nm  : (M, 2) tensor — true emitter positions in nm
    Im0         : emitter intensity in photons/s
    imagej_coords: bool
    """
    device = next(model.parameters()).device
    model.eval()

    true_xy_nm = true_xy_nm.to(device)          # (M, 2)

    # Build SystemConfig to generate one frame
    xy_tensor = true_xy_nm.unsqueeze(0)          # (1, M, 2)
    Im_tensor = torch.full(
        (1, M), float(Im0), dtype=torch.float32, device=device
    )
    emitter = EmitterData(xy=xy_tensor, Im=Im_tensor)
    system  = SystemConfig(
        emitter = emitter,
        psf     = psf,
        camera  = camera,
        noise   = noise,
    )

    with torch.no_grad():
        frame = gauss2d_frame_torch(system)[0]   # (Ky, Kx)

    # Prepare model input
    x_input = frame.flatten().unsqueeze(0)       # (1, Kx*Ky)
    x_input = x_input / (x_input.amax(dim=1, keepdim=True) + 1e-8)

    with torch.no_grad():
        pred_nm = model(x_input)                 # (1, 2*M)

    pred_xy_nm = pred_nm.view(M, 2)              # (M, 2)

    # RMSMD
    D, _ = rmsmd_per_frame(
        pred_xy_nm.unsqueeze(0),                 # (1, M, 2)
        true_xy_nm.unsqueeze(0),                 # (1, M, 2)
    )
    rmsmd_val = D[0].item()

    # Convert nm → pixel coords
    Dx, Dy = camera.Dx, camera.Dy
    true_px = (true_xy_nm / torch.tensor([Dx, Dy], device=device) - 0.5).cpu().numpy()
    pred_px = (pred_xy_nm / torch.tensor([Dx, Dy], device=device) - 0.5).cpu().numpy()
    frame_np = frame.cpu().numpy()

    origin = 'upper' if imagej_coords else 'lower'

    plt.figure(figsize=(5, 5))
    plt.imshow(frame_np, cmap='gray', origin=origin)
    plt.scatter(true_px[:, 0], true_px[:, 1], c='red',   s=40, label='True')
    plt.scatter(pred_px[:, 0], pred_px[:, 1], c='green', s=40, label='Predicted')
    plt.title(f"True vs Predicted Localization\nRMSMD = {rmsmd_val:.2f} nm")
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()


def show8bimage(img0, rescale='Yes', color='gray', addColormap='No'):
    """Display an image in 8-bit style (matches MATLAB show8bimage).

    Parameters
    ----------
    img0        : 2-D array
    rescale     : 'Yes' to stretch to [0,255], else clip
    color       : 'gray' or 'jet'
    addColormap : 'Yes' to append a palette bar
    """
    img0 = np.array(img0, dtype=float)

    if rescale == 'Yes':
        mx, mn = img0.max(), img0.min()
        if mx > mn:
            img1 = (255 * (img0 - mn) / (mx - mn)).astype(np.uint8)
        else:
            img1 = np.zeros_like(img0, dtype=np.uint8)
    else:
        img1 = img0.astype(np.uint8)

    if addColormap == 'Yes':
        img2 = zoom(img1, zoom=(16, 16), order=0)
        R, C = img2.shape
        bar_width = int(np.ceil(C / 20))
        palette = np.floor(255 * np.linspace(1, 0, R)).astype(np.uint8)
        palette = np.tile(palette.reshape(R, 1), (1, bar_width))
        img2 = np.hstack([img2, palette])
    else:
        img2 = img1

    cmap = plt.cm.gray if color == 'gray' else plt.cm.jet
    plt.imshow(img2, cmap=cmap, vmin=0, vmax=255, origin='upper')
    plt.axis('off')
