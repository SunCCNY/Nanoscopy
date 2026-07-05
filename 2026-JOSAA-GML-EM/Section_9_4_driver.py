"""
Section_9_4_driver.py

Simulation driver for Section 9.4 (single frame, I -> infinity, 3D astigmatic).

Verifies Theorem 2 in 3D: as the average emitter intensity I grows with a
single frame N=1, the global maximum-likelihood (GML) estimator of the 3D
astigmatic SMLM model is efficient, its covariance approaching F(I)^{-1}, and
it coincides with the unbiased UGIA-F benchmark. EM-GML is initialized near the
true positions in all three coordinates to locate the GML vector, the analogue
of UGIA-F drawing from the true positions; a position-agnostic initializer such
as SIC (Sun et al., QBI 2014) reaches the GML in practice, so the near-true init
is a stand-in for an achievable one.

Intensity sweep
---------------
The intensity weights beta_m are fixed once, beta_m = draw_m / mean(draw) with
draw_m ~ U[Imin, Imax], so (1/M) sum_m beta_m = 1. The average intensity I is
then swept and the per-emitter intensities are I_m = beta_m I. Unlike Section
9.3, where the N-frame Fisher is N F1 and F1 is computed once, here the
intensity enters F(I) nonlinearly through the per-pixel weight
Qu = Q + (b+G)/I, so the Fisher factor and the CRB are recomputed at each I.

GPU batching (speed)
--------------------
The Monte-Carlo replicates are the batch axis. For each I:
  * the n_mc single frames are generated in ONE astig3d_frame_torch call;
  * EM-GML is ONE astig3d_em_batch_torch call over all n_mc replicates with a
    fixed n_iter (3x3 Newton M-step, all coordinates jittered near truth);
  * UGIA-F is an analytic draw theta0 + N(0, F(I)^{-1}); the per-frame Fisher
    F(I) (KA case, 3M x 3M) is computed once per I, so UGIA-F and the CRB
    reference reuse it.

Metric
------
The estimated and true positions are in known one-to-one correspondence, so the
RMSE is the direct per-emitter root-mean-square error over the three
coordinates, sqrt( mean_{r,m} || xyz[r,m] - xyz0[m] ||^2 ). The single-frame
CRB reference is sqrt( tr(F(I)^{-1}) / M ).

Geometry
--------
FOV [0,Lx]x[0,Ly]x[-Lz,Lz]; emitters in the central block
[200,2200]x[200,2200] nm laterally (margin_px=2 of a 24x24, 100 nm camera) and
[-Lz,Lz]=[-400,400] nm axially. M=160 in 4 um^2 laterally, delta=75 nm in 3D.

Yi Sun / SMLM_Lib
"""

import sys
import os
import math
import time

import torch
import matplotlib.pyplot as plt
import numpy as np

ROOT = r"K:\Research_KB"   # <- change to your actual path
sys.path.append(ROOT)

from SMLM_Lib import (
    Camera, Astigmatic3DPSF, Noise, EmitterData, SystemConfig,
    astig3d_frame_torch, astig3d_em_batch_torch, astig3d_fisher_torch,
    astig3d_snr_torch, sample_delta_separated_3d,
    background_cloud, print_device_info,
)


# ======================================================================
# Fixed configuration (Section 9.4)
# ======================================================================

def build_config_3d(M=40, delta=250.0, K=24, margin_px=2, Lz=400.0, corr_px=10.0,
                    Imin=250000.0, Imax=350000.0,
                    b_range=(0.2, 0.3), G_range=(0.1, 0.2),
                    seed=0, device="cpu", dtype=torch.float32,
                    delta_xy=None, delta_z=None):
    """Build the fixed Section 9.4 configuration (reproducible from seed).

    Returns xyz0 (M,3), beta (M,) with mean 1, the astigmatic PSF, camera, and
    noise. The per-emitter intensities at average intensity I are I_m = beta_m I.
    """
    camera = Camera(Kx=K, Ky=K, Dx=100.0, Dy=100.0, Dt=0.01)
    psf = Astigmatic3DPSF(c=205.0, d=290.0,
                          sigmax0=140.0, Ax=0.05, Bx=0.03,
                          sigmay0=135.0, Ay=-0.01, By=0.02)

    xyz0 = sample_delta_separated_3d(M, camera, delta, Lz, margin_px=margin_px,
                                     seed=seed, device=device,
                                     delta_xy=delta_xy, delta_z=delta_z).to(dtype)

    # Fixed intensity weights beta_m, normalized to mean 1.
    gI = torch.Generator().manual_seed(seed + 1)
    draw = Imin + (Imax - Imin) * torch.rand(M, generator=gI)
    beta = (draw / draw.mean()).to(device, dtype)            # (M,), mean 1

    b_map = background_cloud(camera, b_range[0], b_range[1], corr_px=corr_px,
                             seed=seed + 2, device=device, dtype=dtype)
    gG = torch.Generator().manual_seed(seed + 3)
    G_map = (G_range[0] + (G_range[1] - G_range[0])
             * torch.rand(camera.Ky, camera.Kx, generator=gG)).to(device, dtype)
    noise = Noise(b=b_map, mu=G_map.clone(), G=G_map)         # mu = G
    return xyz0, beta, psf, camera, noise


# ======================================================================
# Fisher-based quantities at intensity I (recomputed per I)
# ======================================================================

def fisher_factor_3d(xyz0, Im, psf, camera, noise):
    """Single-frame 3D position Fisher (KA) at intensities Im.

    Returns L with L @ L^T = F(I)^{-1} (the 3M x 3M position covariance) in the
    interleaved (x0,y0,z0,x1,y1,z1,...) ordering of the FIM, and the per-
    coordinate CRLB variances diag = diag(F(I)^{-1}) in that same order. The
    UGIA-F draw is theta0 + L z, z ~ N(0, I_{3M}); the total, lateral, and axial
    CRB on the RMSE follow from diag as sqrt(sum(.)/M). Computed in float64.
    """
    nb = Noise(noise.b.double(), noise.mu.double(), noise.G.double())
    sys1 = SystemConfig(
        EmitterData(xyz=xyz0.double().unsqueeze(0), Im=Im.double().unsqueeze(0)),
        psf, camera, nb)
    F = astig3d_fisher_torch(sys1, 'KA')                     # (3M, 3M), float64
    crlb = torch.linalg.pinv(F)
    diag = crlb.diagonal().clone()                           # (3M,) per-coord var
    lam, V = torch.linalg.eigh(crlb)
    lam = lam.clamp(min=0.0)
    L = V * lam.sqrt()                                       # L @ L^T = crlb
    return L, diag, lam, V


# ======================================================================
# Near-true (oracle) init, all three coordinates
# ======================================================================

def jitter_init_3d(xyz0, n_mc, near_half=20.0, near_half_xy=None, near_half_z=None):
    """Near-true init with optional anisotropic jitter. The lateral coordinates
    get +-near_half_xy nm and the axial coordinate +-near_half_z nm, each
    defaulting to near_half. The near-true init is what makes EM reach the
    GLOBAL maximum (GML). Because the astigmatic axial mode converges slowly, a
    smaller near_half_z lowers the fixed-iteration axial floor. To keep each
    emitter inside its own basin, the lateral jitter must satisfy
    near_half_xy < delta_xy/(2 sqrt2) and the axial near_half_z < delta_z/2.
    """
    M = xyz0.shape[0]
    device, dtype = xyz0.device, xyz0.dtype
    hxy = near_half if near_half_xy is None else near_half_xy
    hz = near_half if near_half_z is None else near_half_z
    scale = torch.tensor([hxy, hxy, hz], device=device, dtype=dtype)
    off = scale * (2.0 * torch.rand(n_mc, M, 3, device=device, dtype=dtype) - 1.0)
    return xyz0.unsqueeze(0) + off                            # (n_mc, M, 3)


# ======================================================================
# RMSE (known one-to-one correspondence)
# ======================================================================

def rmse_split(xyz0, xyz_batch):
    """RMSE with KNOWN correspondence, returned as (total, lateral, axial) nm:
        total   = sqrt( mean_{r,m} (dx^2 + dy^2 + dz^2) ),
        lateral = sqrt( mean_{r,m} (dx^2 + dy^2) ),
        axial   = sqrt( mean_{r,m} dz^2 ),
    so total^2 = lateral^2 + axial^2. The split exposes the slow axial
    convergence behind the high-I EM floor."""
    d = xyz_batch - xyz0.unsqueeze(0)                         # (n_mc, M, 3)
    e_xy = d[..., 0] ** 2 + d[..., 1] ** 2
    e_z = d[..., 2] ** 2
    rmse_xy = e_xy.mean().sqrt().item()
    rmse_z = e_z.mean().sqrt().item()
    rmse_all = (e_xy + e_z).mean().sqrt().item()
    return rmse_all, rmse_xy, rmse_z


# ======================================================================
# Main experiment
# ======================================================================

def bias_known_3d(xyz0, xyz_batch):
    """RMS bias over emitters, sqrt( mean_m || mean_r xyz[r,m] - xyz0[m] ||^2 ).

    The Monte-Carlo mean estimate minus the truth, in nm, in 3D. Property 1.
    """
    mean_est = xyz_batch.mean(0)                             # (M, 3)
    err2 = ((mean_est - xyz0) ** 2).sum(-1)                  # (M,)
    return err2.mean().sqrt().item()


def run_section_9_4(M=40, delta=250.0, K=24, margin_px=2, Lz=400.0, corr_px=10.0,
                    I_list=(1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8, 3e8, 1e9, 3e9, 1e10),
                    n_mc=20, n_iter=4096, near_half=5.0, GN='N', seed=0,
                    em_dtype=torch.float64, device=None, out_dir="Sec9_4",
                    delta_xy=None, delta_z=None,
                    near_half_xy=None, near_half_z=None):
    """Run the Section 9.4 experiment (GPU-batched) and save the three figures."""
    print_device_info()
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    os.makedirs(out_dir, exist_ok=True)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    xyz0, beta, psf, camera, noise = build_config_3d(
        M=M, delta=delta, K=K, margin_px=margin_px, Lz=Lz, corr_px=corr_px,
        seed=seed, device=device, delta_xy=delta_xy, delta_z=delta_z)

    blk_x = camera.Kx - 2 * margin_px
    blk_y = camera.Ky - 2 * margin_px
    occ_um2 = (blk_x * camera.Dx) * (blk_y * camera.Dy) / 1e6
    print(f"Config: M={M}, delta={delta} nm (3D), frame {camera.Kx}x{camera.Ky} px, "
          f"block {blk_x}x{blk_y} px ({occ_um2:.1f} um^2), z in [-{Lz:g},{Lz:g}] nm, "
          f"lateral density={M / occ_um2:.1f}/um^2")
    hxy = near_half if near_half_xy is None else near_half_xy
    hz = near_half if near_half_z is None else near_half_z
    dxy = delta if delta_xy is None else delta_xy
    dz = delta if delta_z is None else delta_z
    print(f"        n_mc={n_mc}, n_iter={n_iter}, GML init xy +-{hxy:g} nm, "
          f"z +-{hz:g} nm, delta_xy={dxy:g} nm, delta_z={dz:g} nm")

    xyz0d = xyz0.double().reshape(-1)                         # (3M,) interleaved
    rec = {"gml": {}, "ugia": {}, "crb": {}, "snr": {}, "bias": {}}
    recon = {}
    recon_Is = [min(I_list, key=lambda I: abs(I - t)) for t in (1e5, 1e8, 1e10)]
    qq_I = max(I_list)                                        # QQ at the efficient end
    res_qq = diag_qq = lam_qq = V_qq = None

    for I in I_list:
        Im = beta * float(I)                                  # (M,)

        # Fisher factor + per-coordinate CRB at this I (KA; I enters nonlinearly).
        L, diag, lam, Veig = fisher_factor_3d(xyz0, Im, psf, camera, noise)
        L = L.to(device)
        d3 = diag.view(M, 3)
        crb = math.sqrt(diag.sum().item() / M)
        crb_xy = math.sqrt(d3[:, :2].sum().item() / M)
        crb_z = math.sqrt(d3[:, 2].sum().item() / M)

        # Per-emitter SNR of the configuration at this I (Eq. 31).
        sysI = SystemConfig(EmitterData(xyz=xyz0.unsqueeze(0), Im=Im.unsqueeze(0)),
                            psf, camera, noise)
        snr = astig3d_snr_torch(sysI)[0]

        # (a) generate n_mc single frames in ONE batched call.
        em_gen = EmitterData(
            xyz=xyz0.unsqueeze(0).expand(n_mc, M, 3).contiguous(),
            Im=Im.unsqueeze(0).expand(n_mc, M).contiguous())
        V = astig3d_frame_torch(SystemConfig(em_gen, psf, camera, noise))

        # (b) batched EM-GML, fixed n_iter, z clamped to the calibrated depth.
        sys_em = SystemConfig(
            EmitterData(xyz=xyz0.unsqueeze(0), Im=Im.unsqueeze(0)),
            psf, camera, noise)
        xyzi = jitter_init_3d(xyz0, n_mc, near_half=near_half,
                              near_half_xy=near_half_xy, near_half_z=near_half_z)
        xyz_gml = astig3d_em_batch_torch(V, xyzi, sys_em, GN=GN, n_iter=n_iter,
                                         dtype=em_dtype, z_clamp=(-Lz, Lz))

        # (c) analytic UGIA-F draws: theta0 + N(0, F(I)^{-1}).
        z = torch.randn(n_mc, 3 * M, device=device, dtype=torch.float64)
        W = z @ L.t()
        xyzF = (xyz0d + W).reshape(n_mc, M, 3).to(xyz0.dtype)

        g_all, g_xy, g_z = rmse_split(xyz0, xyz_gml)
        u_all, u_xy, u_z = rmse_split(xyz0, xyzF)
        rec["gml"][I] = (g_all, g_xy, g_z)
        rec["ugia"][I] = (u_all, u_xy, u_z)
        rec["crb"][I] = (crb, crb_xy, crb_z)
        rec["snr"][I] = snr.mean().item()
        rec["bias"][I] = bias_known_3d(xyz0, xyz_gml)
        if I in recon_Is:
            recon[I] = {"gml": xyz_gml[0], "ugia": xyzF[0], "frame": V[0]}
        if I == qq_I:
            res_qq = (xyz_gml - xyz0).detach()
            diag_qq = diag.clone(); lam_qq = lam.clone(); V_qq = Veig.clone()
        print(f"  I={I:10.0f}  GML {g_all:7.3f} (xy {g_xy:6.3f} z {g_z:6.3f})  "
              f"UGIA-F {u_all:8.3f}  CRB {crb:8.3f} (xy {crb_xy:6.3f} z {crb_z:6.3f}) "
              f"nm  SNR {rec['snr'][I]:5.1f} dB")

    _make_figures_3d(xyz0, camera, Lz, recon, rec, list(I_list), recon_Is, out_dir,
                     res_qq, qq_I, diag_qq, lam_qq, V_qq, M, n_mc)
    return rec


# ======================================================================
# Figures
# ======================================================================

def _make_figures_3d(xyz0, camera, Lz, recon, rec, I_list, recon_Is, out_dir,
                     res_qq, qq_I, diag_qq, lam_qq, V_qq, M, n_mc):
    t0 = xyz0.detach().cpu().numpy()
    ext = [0, camera.Lx, 0, camera.Ly]
    Is = list(I_list)

    def _save(name, dpi=130):
        # Save each figure as both .png and .svg (vector, for journals).
        base = os.path.join(out_dir, name)
        plt.savefig(base + ".png", dpi=dpi)
        plt.savefig(base + ".svg")

    def _tag(ax, txt, dark=True):
        c, bg = ("w", "k") if dark else ("k", "w")
        ax.text(0.04, 0.95, txt, transform=ax.transAxes, color=c, fontsize=11,
                va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.15", fc=bg, ec="none", alpha=0.5))

    def add_snr_axis(a):
        # Top axis labelling each I with the configuration's mean SNR (dB) at the
        # same log-x positions; a twin axis is used because secondary_xaxis
        # mis-places ticks under a custom map on a log-scaled parent.
        sec = a.twiny()
        sec.set_xscale("log")
        sec.set_xlim(a.get_xlim())
        sec.set_xticks(Is)
        sec.set_xticklabels([f"{rec['snr'][I]:.0f}" for I in Is], fontsize=8)
        sec.minorticks_off()
        sec.set_xlabel("Mean SNR (dB)")

    # Figure (i): 3x3 panel. Rows are I=1e5, 1e8, 1e10; columns are the data
    # frame (x,y), the x-y reconstruction, and the x-z reconstruction. The true
    # positions are the red dots, drawn on top, in every panel.
    rIs = sorted(recon_Is)
    dot_s = 4
    fig, axes = plt.subplots(3, 3, figsize=(11, 10.5))
    for i, I in enumerate(rIs):
        r = recon[I]
        frame = r["frame"].detach().cpu().numpy()
        g = r["gml"].detach().cpu().numpy()
        u = r["ugia"].detach().cpu().numpy()
        ax = axes[i, 0]
        ax.imshow(frame, origin="lower", extent=ext, cmap="viridis")
        ax.scatter(t0[:, 0], t0[:, 1], s=dot_s, c="red", zorder=5)
        ax.set_xlim(0, camera.Lx); ax.set_ylim(0, camera.Ly)
        _tag(ax, f"$I=10^{{{round(math.log10(I))}}}$")
        ax = axes[i, 1]
        ax.scatter(u[:, 0], u[:, 1], s=16, marker="+", c="tab:green", label="UGIA-F")
        ax.scatter(g[:, 0], g[:, 1], s=16, marker="x", c="tab:blue", label="EM-GML")
        ax.scatter(t0[:, 0], t0[:, 1], s=dot_s, c="red", zorder=5, label="True positions")
        ax.set_xlim(0, camera.Lx); ax.set_ylim(0, camera.Ly)
        ax = axes[i, 2]
        ax.scatter(u[:, 0], u[:, 2], s=16, marker="+", c="tab:green")
        ax.scatter(g[:, 0], g[:, 2], s=16, marker="x", c="tab:blue")
        ax.scatter(t0[:, 0], t0[:, 2], s=dot_s, c="red", zorder=5)
        ax.set_xlim(0, camera.Lx); ax.set_ylim(-Lz, Lz)
    _tag(axes[0, 1], "x-y", dark=False)
    _tag(axes[0, 2], "x-z", dark=False)
    for ax in axes.ravel():
        ax.set_xticks([]); ax.set_yticks([])
    # 500 nm lateral scale bar in the I=1e5 frame (top-left), butt cap so it is
    # exactly 5 pixels long.
    bar, xb, yb = 500.0, camera.Lx * 0.55 + 500, camera.Ly * 0.08 - 100
    axes[0, 0].plot([xb, xb + bar], [yb, yb], "-", color="w", lw=3, zorder=6,
                    solid_capstyle="butt")
    axes[0, 0].text(xb + bar / 2, yb + camera.Ly * 0.02, "500 nm", color="w",
                    ha="center", va="bottom", fontsize=8, zorder=6)
    axes[0, 1].legend(loc="lower left", fontsize=7, framealpha=0.6, markerscale=1)
    plt.tight_layout()
    _save("fig_3_frames_and_recon_3d", dpi=140)
    plt.close()

    # ----- Figure 4: total RMSE, bias, and QQ as a 3x1 stack -----
    # Combined for the journal figure budget, mirroring Figure 2 of Section 9.3.
    # (a) total RMSE (x, y, z) versus I for EM-GML and UGIA-F with the CRB; the
    # lateral and axial RMSE are nearly identical, so only the total is shown and
    # the caption notes the equality. (b) bias versus I with the Monte-Carlo
    # floor RMSE/sqrt(L), the measured bias riding the floor so the RMSE is
    # variance-dominated. (c) whitened-residual QQ at qq_I, the residual being
    # z = F(I)^{1/2}(xyz_hat - xyz0) ~ N(0, I) under Theorem 2, verifying
    # normality and per-coordinate efficiency.
    print(f"  lateral vs axial EM-GML RMSE at I={Is[-1]:.0f}: "
          f"xy={rec['gml'][Is[-1]][1]:.3f} nm, z={rec['gml'][Is[-1]][2]:.3f} nm")
    fig, axes = plt.subplots(3, 1, figsize=(6.0, 14.0))
    ax_rmse, ax_bias, ax_qq = axes

    # (a) total RMSE
    ax_rmse.loglog(Is, [rec["ugia"][I][0] for I in Is], "g+-", label="UGIA-F")
    ax_rmse.loglog(Is, [rec["gml"][I][0] for I in Is], "bx-", label="EM-GML")
    ax_rmse.loglog(Is, [rec["crb"][I][0] for I in Is], "k--", label="CRB")
    ax_rmse.set_xlabel("Average intensity I (photons/s)")
    ax_rmse.set_ylabel("RMSE (nm)")
    ax_rmse.legend(loc="upper right")
    ax_rmse.grid(True, which="both", alpha=0.3)
    add_snr_axis(ax_rmse)

    # (b) bias versus I with the Monte-Carlo floor
    bias = [rec["bias"][I] for I in Is]
    floor = [rec["gml"][I][0] / math.sqrt(n_mc) for I in Is]
    ax_bias.loglog(Is, [rec["gml"][I][0] for I in Is], "bx--", alpha=0.4, label="RMSE")
    ax_bias.loglog(Is, bias, "bo-", label=r"bias $\|\overline{\hat\theta}-\theta_0\|$")
    ax_bias.loglog(Is, floor, "k:", label=r"Monte-Carlo floor RMSE$/\sqrt{L}$")
    ax_bias.set_xlabel("Average intensity I (photons/s)")
    ax_bias.set_ylabel("Error (nm)")
    ax_bias.legend(loc="upper right")
    ax_bias.grid(True, which="both", alpha=0.3)
    add_snr_axis(ax_bias)

    # (c) whitened-residual QQ at qq_I
    if res_qq is not None:
        res = res_qq.reshape(res_qq.shape[0], -1).double()   # (n_mc, 3M), nm
        lam_c = lam_qq.clamp(min=lam_qq.max() * 1e-12)
        Z = (res @ V_qq) / lam_c.sqrt()                      # (n_mc, 3M) ~ N(0,1)
        zf = Z.reshape(-1).detach().cpu().numpy()
        zf.sort()
        n = zf.size
        pp = (np.arange(n) + 0.5) / n
        theo = torch.special.ndtri(torch.from_numpy(pp)).numpy()
        d3 = diag_qq.view(M, 3)                              # (M,3) per-coord var
        std = res_qq.double() / d3.unsqueeze(0).sqrt()       # (n_mc, M, 3)
        ratio_x = (std[..., 0] ** 2).mean().item()
        ratio_y = (std[..., 1] ** 2).mean().item()
        ratio_z = (std[..., 2] ** 2).mean().item()
        print(f"  whitened-residual variance ratio at I={qq_I:.0f}: "
              f"x={ratio_x:.3f}, y={ratio_y:.3f}, z={ratio_z:.3f} (ideal 1.000)")
        ax_qq.plot(theo, zf, ".", ms=2, color="tab:blue")
        lim = [min(theo[0], zf[0]), max(theo[-1], zf[-1])]
        ax_qq.plot(lim, lim, "k--", lw=1, label="Ideal $N(0,1)$")
        ax_qq.set_xlabel("Standard normal quantiles")
        ax_qq.set_ylabel("Whitened residual quantiles")
        ax_qq.legend(loc="lower right")
        ax_qq.grid(True, alpha=0.3)

    # subpanel labels (a)-(c) inside the panels
    for a, lab in ((ax_rmse, "(a)"), (ax_bias, "(b)"), (ax_qq, "(c)")):
        a.text(0.03, 0.5, lab, transform=a.transAxes, fontsize=12,
               fontweight="bold", va="center", ha="left",
               bbox=dict(boxstyle="round,pad=0.15", fc="w", ec="none", alpha=0.7))

    fig.tight_layout()
    _save("fig_4_rmse_bias_qq")
    plt.close()

if __name__ == "__main__":
    torch.set_float32_matmul_precision('high')
    t_start = time.perf_counter()
    run_section_9_4()
    print(f"Total time {time.perf_counter() - t_start:.1f} s")
