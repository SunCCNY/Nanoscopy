"""
Section_9_5_driver.py

Simulation driver for Section 9.5 (accuracy versus emitter density at N = 1).

Sweeps the emitter count M inside a fixed central block, so the density is
M / block_area, and reports ARMSMD-P at a single frame (N = 1) for EM-GML,
UGIA-F, and the Cramer-Rao bound. Reproduces the property of Sun et al.
(QBI 2014): EM-GML accuracy grows roughly linearly with density while UGIA-F
and the CRB grow exponentially, the regime that matters for spatiotemporal
super-resolution.

EM-GML is initialized near the true positions to locate the GML vector, the
analogue of UGIA-F drawing from the true positions; a position-agnostic
initializer such as SIC (Sun et al., QBI 2014) reaches the GML in practice, so
the near-true init is a stand-in for an achievable one.

The configuration and helpers below are identical to Section_9_3_driver.py
(kept in sync by copy), so the two sections use exactly the same system.

Yi Sun / SMLM_Lib
"""

import sys
import os
import math
import time

import torch
import matplotlib.pyplot as plt

ROOT = r"K:\Research_KB"   # <- change to your actual path
sys.path.append(ROOT)

from SMLM_Lib import (
    Camera, Gaussian2DPSF, Noise, EmitterData, SystemConfig,
    gauss2d_frame_torch, gauss2d_em_batch_torch, gauss2d_fisher_torch,
    sample_delta_separated, background_cloud, print_device_info,
)


# ======================================================================
# Fixed configuration (identical to Section 9.3)
# ======================================================================

def build_config(M=160, delta=100.0, corr_px=10.0, K=24, margin_px=2,
                 Imin=250000.0, Imax=350000.0,
                 b_range=(4.0, 6.0), G_range=(2.0, 4.0),
                 seed=0, device="cpu", dtype=torch.float32):
    """Build the fixed configuration (reproducible from seed)."""
    camera = Camera(Kx=K, Ky=K, Dx=100.0, Dy=100.0, Dt=0.01)
    psf = Gaussian2DPSF(sigma=108.81)

    theta0 = sample_delta_separated(M, camera, delta, margin_px=margin_px,
                                    seed=seed, device=device).to(dtype)
    gI = torch.Generator().manual_seed(seed + 1)
    Im = (Imin + (Imax - Imin) * torch.rand(M, generator=gI)).to(device, dtype)

    b_map = background_cloud(camera, b_range[0], b_range[1], corr_px=corr_px,
                             seed=seed + 2, device=device, dtype=dtype)
    gG = torch.Generator().manual_seed(seed + 3)
    G_map = (G_range[0] + (G_range[1] - G_range[0])
             * torch.rand(camera.Ky, camera.Kx, generator=gG)).to(device, dtype)
    noise = Noise(b=b_map, mu=G_map.clone(), G=G_map)   # mu = G
    return theta0, Im, psf, camera, noise


# ======================================================================
# Fisher-based quantities (UGIA-F factor + CRB)
# ======================================================================

def fisher_factor(theta0, Im, psf, camera, noise):
    """Single-frame position Fisher (KA), its sqrt-covariance factor, and trace.

    Returns L1 with L1 @ L1^T = F^{-1} (the per-frame position covariance) in
    the interleaved (x0,y0,x1,y1,...) ordering of the FIM, and tr(F^{-1}). The
    UGIA-F draw at N = 1 is theta0 + L1 @ z, and the ARMSMD-P CRB is
    sqrt(tr(F^{-1})/M). Computed in float64 for stability.
    """
    nb = Noise(noise.b.double(), noise.mu.double(), noise.G.double())
    sys1 = SystemConfig(
        EmitterData(xy=theta0.double().unsqueeze(0), Im=Im.double().unsqueeze(0)),
        psf, camera, nb)
    F = gauss2d_fisher_torch(sys1, 'KA')                 # (2M, 2M), float64
    crlb1 = torch.linalg.pinv(F)
    lam, V = torch.linalg.eigh(crlb1)
    lam = lam.clamp(min=0.0)
    L1 = V * lam.sqrt()                                  # L1 @ L1^T = crlb1
    return L1, crlb1.diagonal().sum().item()


# ======================================================================
# Initialisation + metric
# ======================================================================

def jitter_init(theta0, n_mc, near_half=25.0):
    """Near-true (oracle) init: true positions plus uniform jitter in
    +-near_half nm per coordinate, matching the EM-GML benchmark convention.

    The near-true init is what makes EM reach the GLOBAL maximum (GML). At high
    emitter density the basins of attraction are narrow, so near_half must be
    smaller than half the minimum separation, otherwise an emitter is initialised
    inside a neighbour's basin and the EM converges to a spurious maximum.

    A position-agnostic initializer such as successive interference cancellation
    (SIC; Sun et al., QBI 2014) brings the EM to the GML in practice, so this
    near-true init is a stand-in for an achievable initialization, the EM-GML
    analogue of UGIA-F drawing from the true positions.
    """
    M = theta0.shape[0]
    device, dtype = theta0.device, theta0.dtype
    off = near_half * (2.0 * torch.rand(n_mc, M, 2, device=device, dtype=dtype) - 1.0)
    return theta0.unsqueeze(0) + off                      # (n_mc, M, 2)


def armsmd_p_known(theta0, xy_batch):
    """ARMSMD-P when emitter correspondence is KNOWN (UGIA-F and EM-GML).

    The UGIA-F draw is per emitter and the near-true GML init pins estimate m to
    emitter m, so estimate m corresponds to true emitter m. The RMSMD-P is then
    the direct per-emitter RMSE, averaged over replicates and emitters:
        sqrt( mean_{r,m} || xy[r,m] - theta0[m] ||^2 ).
    This is the exact RMSMD-P for a correct partition and avoids partition_x,
    whose pooled assignment mis-matches at high density.
    """
    err2 = ((xy_batch - theta0.unsqueeze(0)) ** 2).sum(-1)   # (n_mc, M)
    return err2.mean().sqrt().item()


# ======================================================================
# Section 9.5 experiment: accuracy versus density at N = 1
# ======================================================================

def run_section_9_5(M_list=(1, 2, 4, 8, 16, 32, 48, 80, 120, 160, 200),
                    delta=100.0, K=24, margin_px=2, n_mc=20, n_iter=2048,
                    near_half=25.0, GN='N', seed=0, em_dtype=torch.float64,
                    device=None, out_dir="Sec9_5"):
    """ARMSMD-P versus emitter density at N = 1 (single frame).

    Sweeps the emitter count M inside the fixed central block, so the density is
    M / block_area. Shows EM-GML growing roughly linearly with density while
    UGIA-F and the CRB grow exponentially (Sun et al., QBI 2014).
    """
    print_device_info()
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    os.makedirs(out_dir, exist_ok=True)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    dens, gml, ugia, crb = [], [], [], []
    for M in M_list:
        theta0, Im, psf, camera, noise = build_config(
            M=M, delta=delta, K=K, margin_px=margin_px, seed=seed, device=device)
        L1, tr1 = fisher_factor(theta0, Im, psf, camera, noise)
        L1 = L1.to(device)
        theta0d = theta0.double().reshape(-1)
        blk_x = camera.Kx - 2 * margin_px
        blk_y = camera.Ky - 2 * margin_px
        occ_um2 = (blk_x * camera.Dx) * (blk_y * camera.Dy) / 1e6
        D = M / occ_um2

        # single frame (N = 1)
        em_gen = EmitterData(xy=theta0.unsqueeze(0).expand(n_mc, M, 2).contiguous(),
                             Im=Im.unsqueeze(0).expand(n_mc, M).contiguous())
        Vp = gauss2d_frame_torch(SystemConfig(em_gen, psf, camera, noise))
        sys_em = SystemConfig(EmitterData(xy=theta0.unsqueeze(0), Im=Im.unsqueeze(0)),
                              psf, camera, noise)
        xyi = jitter_init(theta0, n_mc, near_half=near_half)
        xy_gml = gauss2d_em_batch_torch(Vp, xyi, sys_em, GN=GN, n_iter=n_iter,
                                        dtype=em_dtype)
        z = torch.randn(n_mc, 2 * M, device=device, dtype=torch.float64)
        xyF = (theta0d + z @ L1.t()).reshape(n_mc, M, 2).to(theta0.dtype)

        dens.append(D)
        gml.append(armsmd_p_known(theta0, xy_gml))
        ugia.append(armsmd_p_known(theta0, xyF))
        crb.append(math.sqrt(tr1 / M))
        print(f"  M={M:4d}  density={D:5.1f}/um^2  GML={gml[-1]:8.2f}  "
              f"UGIA-F={ugia[-1]:8.2f}  CRB={crb[-1]:8.2f}  nm")

    _make_density_figure(dens, gml, ugia, crb, out_dir)
    return {"density": dens, "gml": gml, "ugia": ugia, "crb": crb}


def _make_density_figure(dens, gml, ugia, crb, out_dir):
    import numpy as np
    D = np.array(dens); g = np.array(gml); u = np.array(ugia); c = np.array(crb)
    ag, bg = np.polyfit(D, g, 1)                      # EM-GML: linear in D
    sc, ic = np.polyfit(D, np.log(c), 1)             # CRB: log-linear (exponential)
    Ac, rc = math.exp(ic), math.exp(sc)
    Dg = np.linspace(D.min(), D.max(), 100)

    def _save(name, dpi=130):
        # Save the figure as both .png and .svg (vector, for journals).
        base = os.path.join(out_dir, name)
        plt.savefig(base + ".png", dpi=dpi)
        plt.savefig(base + ".svg")

    plt.figure(figsize=(6.6, 4.9))
    plt.semilogy(D, u, "g+", ms=10, label="UGIA-F")
    plt.semilogy(D, c, "k^", ms=6, mfc="none", label="CRB")
    plt.semilogy(Dg, Ac * rc ** Dg, "k--", lw=1,
                 label=rf"CRB fit {Ac:.2f}$\times${rc:.3f}$^D$")
    plt.semilogy(D, g, "bo", ms=6, label="EM-GML")
    plt.semilogy(Dg, ag * Dg + bg, "b-", lw=1,
                 label=rf"GML fit {ag:.2f}$\times$D+{bg:.2f}")
    plt.xlabel(r"Emitter density $D$ (emitters / $\mu$m$^2$)")
    plt.ylabel("RMSE (nm)")
    plt.legend(fontsize=8)
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    _save("fig_5_density_sweep")
    plt.close()


if __name__ == "__main__":
    torch.set_float32_matmul_precision('high')
    t_start = time.perf_counter()
    run_section_9_5()
    print(f"Total time {time.perf_counter() - t_start:.1f} s")
