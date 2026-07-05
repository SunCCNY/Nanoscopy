"""
Section_9_3_driver.py

Simulation driver for Section 9.3 (multiple frames, N -> infinity).

Verifies Theorem 1: the N-frame global maximum-likelihood (GML) estimator is
consistent and efficient (covariance (N F)^{-1}), coinciding with the UGIA-F
benchmark at low density and, being biased, falling below the unbiased CRB at
high density. EM-GML is initialized near the true positions to locate the GML
vector, the analogue of UGIA-F drawing from the true positions; a position-
agnostic initializer such as SIC (Sun et al., QBI 2014) reaches the GML in
practice, so the near-true init is a stand-in for an achievable one.

Sum-of-frames approach
----------------------
The N i.i.d. frames share the same mean v(k), so their sum
    V_+(k) = sum_n V^(n)(k)  ~  Poisson(N v(k))
is a sufficient statistic. The N-frame GML is one EM run on V_+ with the signal
and noise means scaled by N (Im, b, mu, G all x N), so the Poisson mean is
N v(k). Because the N frames are i.i.d., V_+ is sampled directly in one batched
gauss2d_frame_torch call on the N-scaled system rather than summing N frames.

GPU batching (speed)
--------------------
The Monte-Carlo replicates are the batch axis. For each N:
  * the n_mc summed frames are generated in ONE gauss2d_frame_torch call;
  * EM-GML is ONE gauss2d_em_batch_torch call over all n_mc replicates with a
    fixed n_iter (no per-frame convergence loop);
  * UGIA-F is an analytic draw theta0 + N(0, (N F)^{-1}); the per-frame Fisher
    F is computed once (it depends only on the fixed configuration), so UGIA-F
    and the CRB reference both reuse it. This also supports the unequal-
    intensity ('KA') case, which the batched UGIA-F routine does not.

Noise model
-----------
Background b_k: a smooth autofluorescence cloud (one static map). Readout G_k:
per-pixel independent (sCMOS), mean mu_k = G_k. Both maps fixed across frames.

Yi Sun / SMLM_Lib
"""

import sys
import os
import math
import time

import numpy as np
import torch
import matplotlib.pyplot as plt

ROOT = r"K:\Research_KB"   # <- change to your actual path
sys.path.append(ROOT)

from SMLM_Lib import (
    Camera, Gaussian2DPSF, Noise, EmitterData, SystemConfig,
    gauss2d_frame_torch, gauss2d_em_batch_torch, gauss2d_fisher_torch,
    gauss2d_snr_torch, sample_delta_separated,
    background_cloud, print_device_info,
)


# ======================================================================
# Fixed configuration (Section 9.3)
# ======================================================================

def build_config(M=80, delta=120.0, corr_px=10.0, K=24, margin_px=2,
                 Imin=250000.0, Imax=350000.0,
                 b_range=(4.0, 6.0), G_range=(2.0, 4.0),
                 seed=0, device="cpu", dtype=torch.float32):
    """Build the fixed Section 9.3 configuration (reproducible from seed)."""
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


def scaled_noise(noise, s):
    """Noise with all densities scaled by s (the N-frame sum)."""
    return Noise(b=noise.b * s, mu=noise.mu * s, G=noise.G * s)


# ======================================================================
# Fisher-based quantities (computed once: UGIA-F factor + CRB)
# ======================================================================

def fisher_factor(theta0, Im, psf, camera, noise):
    """Single-frame position Fisher (KA), its sqrt-covariance factor, and trace.

    Returns L1 with L1 @ L1^T = F^{-1} (the per-frame position covariance) in
    the interleaved (x0,y0,x1,y1,...) ordering of the FIM, and tr(F^{-1}). The
    N-frame covariance is F^{-1}/N, so the UGIA-F factor is L1/sqrt(N) and the
    ARMSMD-P CRB is sqrt(tr(F^{-1})/(N M)). Computed in float64 for stability.
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
    return L1, crlb1.diagonal().sum().item(), crlb1, lam, V


# ======================================================================
# Initialisations
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


# ======================================================================
# ARMSMD-P
# ======================================================================

def armsmd_p_known(theta0, xy_batch):
    """ARMSMD-P when emitter correspondence is KNOWN (UGIA-F and EM-GML).

    The UGIA-F draw is per emitter and the near-true GML init pins estimate m to
    emitter m, so estimate m corresponds to true emitter m. The RMSMD-P is then
    the direct per-emitter RMSE, averaged over replicates and emitters:
        sqrt( mean_{r,m} || xy[r,m] - theta0[m] ||^2 ).
    This is the exact RMSMD-P for a correct partition and avoids partition_x,
    whose pooled assignment mis-matches at high density (it inflated UGIA-F well
    above the CRB on GPU even though UGIA-F is, by construction, on the bound).
    """
    err2 = ((xy_batch - theta0.unsqueeze(0)) ** 2).sum(-1)   # (n_mc, M)
    return err2.mean().sqrt().item()


def bias_known(theta0, xy_batch):
    """RMS bias over emitters, sqrt( mean_m || mean_r xy[r,m] - theta0[m] ||^2 ).

    The Monte-Carlo mean estimate minus the truth, in nm. Under asymptotic
    unbiasedness this falls as O(1/N), faster than the O(1/sqrt N) standard
    deviation, so the RMSE is variance-dominated. Property 1 of Section 9.1.
    """
    mean_est = xy_batch.mean(0)                              # (M, 2)
    err2 = ((mean_est - theta0) ** 2).sum(-1)                # (M,)
    return err2.mean().sqrt().item()



# ======================================================================
# Main experiment
# ======================================================================

def run_section_9_3(M=80, delta=120.0, corr_px=10.0, K=24, margin_px=2,
                    N_list=(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000),
                    n_mc=20, n_iter=2048, near_half=5.0, GN='N', seed=0,
                    em_dtype=torch.float64, device=None, out_dir="Sec9_3",
                    qq_N=None):
    """Run the Section 9.3 experiment (GPU-batched) and save the three figures."""
    print_device_info()
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    os.makedirs(out_dir, exist_ok=True)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    theta0, Im, psf, camera, noise = build_config(
        M=M, delta=delta, corr_px=corr_px, K=K, margin_px=margin_px,
        seed=seed, device=device)

    # Information-sufficient SNR of the configuration (Eq. 29).
    sys0 = SystemConfig(EmitterData(xy=theta0.unsqueeze(0), Im=Im.unsqueeze(0)),
                        psf, camera, noise)
    snr = gauss2d_snr_torch(sys0)[0]
    blk_x = camera.Kx - 2 * margin_px
    blk_y = camera.Ky - 2 * margin_px
    occ_um2 = (blk_x * camera.Dx) * (blk_y * camera.Dy) / 1e6
    print(f"Config: M={M}, delta={delta} nm, frame {camera.Kx}x{camera.Ky} px, "
          f"emitter block {blk_x}x{blk_y} px ({occ_um2:.1f} um^2), "
          f"density={M / occ_um2:.1f}/um^2")
    print(f"        n_mc={n_mc}, n_iter={n_iter}, GML init +-{near_half:g} nm")
    print(f"Per-emitter SNR (Eq.29): mean {snr.mean():.2f} dB, "
          f"range [{snr.min():.2f}, {snr.max():.2f}] dB")

    # Fisher-based UGIA-F factor + CRB trace + per-coord CRLB (computed once).
    L1, tr1, crlb1, lam, Veig = fisher_factor(theta0, Im, psf, camera, noise)
    L1 = L1.to(device)
    crlb1, lam, Veig = crlb1.to(device), lam.to(device), Veig.to(device)
    theta0d = theta0.double().reshape(-1)                 # (2M,) interleaved

    if qq_N is None:
        qq_N = max(N_list)                                # whitened-residual figure at largest N
    res_qq = None

    rec = {"gml": {}, "ugia": {}, "crb": {}, "bias": {}}
    recon = {}

    for N in N_list:
        crb = math.sqrt(tr1 / (N * M))
        Im_N = N * Im
        noise_N = scaled_noise(noise, N)

        # (a) generate n_mc summed frames in ONE batched call.
        em_gen = EmitterData(xy=theta0.unsqueeze(0).expand(n_mc, M, 2).contiguous(),
                             Im=Im_N.unsqueeze(0).expand(n_mc, M).contiguous())
        Vp = gauss2d_frame_torch(SystemConfig(em_gen, psf, camera, noise_N))

        # (b) batched EM-GML, fixed n_iter.
        sys_em = SystemConfig(EmitterData(xy=theta0.unsqueeze(0), Im=Im_N.unsqueeze(0)),
                              psf, camera, noise_N)
        xyi_gml = jitter_init(theta0, n_mc, near_half=near_half)
        xy_gml = gauss2d_em_batch_torch(Vp, xyi_gml, sys_em, GN=GN, n_iter=n_iter,
                                        dtype=em_dtype)

        # (c) analytic UGIA-F draws: theta0 + N(0, (N F)^{-1}).
        z = torch.randn(n_mc, 2 * M, device=device, dtype=torch.float64)
        W = (z @ L1.t()) / math.sqrt(N)
        xyF = (theta0d + W).reshape(n_mc, M, 2).to(theta0.dtype)

        rec["gml"][N] = armsmd_p_known(theta0, xy_gml)
        rec["ugia"][N] = armsmd_p_known(theta0, xyF)
        rec["crb"][N] = crb
        rec["bias"][N] = bias_known(theta0, xy_gml)
        if N == qq_N:
            res_qq = (xy_gml.double() - theta0.double())  # (n_mc, M, 2), nm
        if N in (1, 2, 5, 10, 100, 1000):
            recon[N] = {"gml": xy_gml[0], "ugia": xyF[0], "frame": Vp[0]}
        print(f"  N={N:5d}  GML={rec['gml'][N]:7.3f}  "
              f"UGIA-F={rec['ugia'][N]:7.3f}  CRB={crb:7.3f}  nm")

    # pure-noise frame (background + readout, no emitter signal) for figure (i).
    em_zero = EmitterData(xy=theta0.unsqueeze(0),
                          Im=torch.zeros_like(Im).unsqueeze(0))
    recon["noise"] = gauss2d_frame_torch(
        SystemConfig(em_zero, psf, camera, noise))[0]

    _make_figures(theta0, camera, recon, rec, N_list, out_dir,
                  res_qq, qq_N, crlb1, lam, Veig, M, n_mc)
    return rec


# ======================================================================
# Figures
# ======================================================================

def _make_figures(theta0, camera, recon, rec, N_list, out_dir,
                  res_qq, qq_N, crlb1, lam, Veig, M, n_mc):
    t0 = theta0.detach().cpu().numpy()
    ext = [0, camera.Lx, 0, camera.Ly]

    def _save(name, dpi=130):
        # Save each figure as both .png and .svg (vector, for journals).
        base = os.path.join(out_dir, name)
        plt.savefig(base + ".png", dpi=dpi)
        plt.savefig(base + ".svg")

    # Figure (i): combined 2x4 panel. Row 1 is the pure-noise frame and the
    # summed frames at N=1,10,100; row 2 is the summed frame at N=1000 and the
    # reconstructions at N=1,10,1000. Each frame uses ordinary per-panel
    # autoscaling, so its brightest pixel maps to the top of the colormap; the
    # noise still visibly shrinks as N grows because the summed frame's relative
    # noise falls as 1/sqrt(N). The standalone true-position map is dropped, as
    # every panel already carries the true positions.
    frame_Ns = [1, 10, 100, 1000]
    if "noise" in recon and all(N in recon for N in frame_Ns):
        frames = {N: recon[N]["frame"].detach().cpu().numpy() for N in frame_Ns}
        noise_img = recon["noise"].detach().cpu().numpy()
        imkw = dict(origin="lower", extent=ext, cmap="viridis")

        def _tag(ax, txt, dark=True):
            c, bg = ("w", "k") if dark else ("k", "w")
            ax.text(0.04, 0.95, txt, transform=ax.transAxes, color=c,
                    fontsize=11, va="top", ha="left",
                    bbox=dict(boxstyle="round,pad=0.15", fc=bg, ec="none",
                              alpha=0.5))

        dot_s = 3                          # true-position red dot, same in all panels

        fig, axes = plt.subplots(2, 4, figsize=(12, 6.2))
        for ax in axes.ravel():
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_xlim(0, camera.Lx); ax.set_ylim(0, camera.Ly)
            ax.set_aspect("equal")

        # Row 1: pure-noise frame, then summed frames N=1, 10, 100. The same
        # small red dot marks the true positions in every panel, drawn on top.
        row1 = [("noise", noise_img), ("N=1", frames[1]),
                ("N=10", frames[10]), ("N=100", frames[100])]
        for ax, (txt, img) in zip(axes[0], row1):
            ax.imshow(img, **imkw)
            ax.scatter(t0[:, 0], t0[:, 1], s=dot_s, c="red", zorder=5,
                       label="True positions")
            _tag(ax, txt)
        # Legend markers at the same size as in the panels (markerscale=1), and
        # a 500 nm scale bar, which the dropped axis ticks would otherwise give.
        axes[0, 0].legend(loc="lower left", fontsize=8, framealpha=0.6,
                          markerscale=1)
        bar, xb, yb = 500.0, camera.Lx * 0.60 + 300, camera.Ly * 0.07 - 100
        axes[0, 0].plot([xb, xb + bar], [yb, yb], "-", color="w", lw=3, zorder=6,
                        solid_capstyle="butt")
        axes[0, 0].text(xb + bar / 2, yb + camera.Ly * 0.02, "500 nm",
                        color="w", ha="center", va="bottom", fontsize=8, zorder=6)

        # Row 2: summed frame N=1000, then reconstructions N=1, 10, 1000. The
        # true red dot is drawn on top of the estimate markers (higher zorder).
        axes[1, 0].imshow(frames[1000], **imkw)
        axes[1, 0].scatter(t0[:, 0], t0[:, 1], s=dot_s, c="red", zorder=5)
        _tag(axes[1, 0], "N=1000")
        for ax, N in zip(axes[1, 1:], [1, 10, 1000]):
            ax.scatter(t0[:, 0], t0[:, 1], s=dot_s, c="red", zorder=5,
                       label="True positions")
            for key, mk, col, lab in (("ugia", "+", "tab:green", "UGIA-F"),
                                      ("gml", "x", "tab:blue", "EM-GML")):
                e = recon[N][key].detach().cpu().numpy()
                ax.scatter(e[:, 0], e[:, 1], s=16, marker=mk, c=col, label=lab)
            _tag(ax, f"N={N}", dark=False)
        axes[1, 1].legend(loc="lower left", fontsize=7, framealpha=0.6,
                          markerscale=1)

        plt.tight_layout()
        plt.subplots_adjust(wspace=0.04, hspace=0.04)
        _save("fig_1_frames_and_recon", dpi=140)
        plt.close()

    # ----- Figure 2: RMSE, bias, and whitened-residual QQ as a 3x1 stack -----
    # Combined for the journal figure budget. (a) RMSE versus N verifies
    # consistency and efficiency, UGIA-F, EM-GML, and CRB meeting as N grows;
    # (b) bias versus N with the Monte-Carlo floor RMSE/sqrt(L) shows the bias
    # rides the floor, so the RMSE is variance-dominated, property 1;
    # (c) whitened-residual QQ at qq_N verifies normality, property 3, and
    # per-coordinate efficiency, property 2, the whitened residual being
    # z = (N F)^{1/2}(theta_hat - theta0) ~ N(0, I) under Theorem 1.
    Ns = list(N_list)
    fig, axes = plt.subplots(3, 1, figsize=(6.0, 14.0))

    # (a) RMSE versus N
    ax = axes[0]
    ax.loglog(Ns, [rec["ugia"][N] for N in Ns], "g+-", label="UGIA-F")
    ax.loglog(Ns, [rec["gml"][N] for N in Ns], "bx-", label="EM-GML")
    ax.loglog(Ns, [rec["crb"][N] for N in Ns], "k--",
              label=r"CRB $\sqrt{\mathrm{tr}(F^{-1})/(NM)}$")
    ax.set_xlabel("Number of frames N")
    ax.set_ylabel("RMSE (nm)")
    ax.legend(loc="upper right")
    ax.grid(True, which="both", alpha=0.3)

    # (b) bias versus N with the Monte-Carlo floor RMSE/sqrt(L)
    bias = [rec["bias"][N] for N in Ns]
    floor = [rec["gml"][N] / math.sqrt(n_mc) for N in Ns]
    ax = axes[1]
    ax.loglog(Ns, [rec["gml"][N] for N in Ns], "bx--", alpha=0.4, label="RMSE")
    ax.loglog(Ns, bias, "bo-", label=r"bias $\|\overline{\hat\theta}-\theta_0\|$")
    ax.loglog(Ns, floor, "k:", label=r"Monte-Carlo floor RMSE$/\sqrt{L}$")
    ax.set_xlabel("Number of frames N")
    ax.set_ylabel("Error (nm)")
    ax.legend(loc="upper right")
    ax.grid(True, which="both", alpha=0.3)

    # (c) whitened-residual QQ at qq_N
    ax = axes[2]
    if res_qq is not None:
        res = res_qq.reshape(res_qq.shape[0], -1)            # (n_mc, 2M), nm
        lam_c = lam.clamp(min=lam.max() * 1e-12)
        Z = (res @ Veig) / lam_c.sqrt() * math.sqrt(qq_N)    # (n_mc, 2M)
        z = Z.reshape(-1).detach().cpu().numpy()
        z.sort()
        n = z.size
        pp = (np.arange(n) + 0.5) / n
        theo = torch.special.ndtri(torch.from_numpy(pp)).numpy()
        diag = crlb1.diagonal().reshape(M, 2)                # per-frame, (M,2)
        std = res_qq / (diag.unsqueeze(0) / qq_N).sqrt()     # (n_mc, M, 2)
        ratio_x = (std[..., 0] ** 2).mean().item()
        ratio_y = (std[..., 1] ** 2).mean().item()
        ax.plot(theo, z, ".", ms=2, color="tab:blue")
        lim = [min(theo[0], z[0]), max(theo[-1], z[-1])]
        ax.plot(lim, lim, "k--", lw=1, label="Ideal $N(0,1)$")
        ax.set_xlabel("Standard normal quantiles")
        ax.set_ylabel("Whitened residual quantiles")
        ax.legend(loc="lower right")
        ax.grid(True, alpha=0.3)
        print(f"  whitened-residual variance ratio at N={qq_N}: "
              f"x={ratio_x:.3f}, y={ratio_y:.3f} (ideal 1.000)")

    # subpanel labels (a)-(c) inside the panels
    for a, lab in zip(axes, ["(a)", "(b)", "(c)"]):
        a.text(0.03, 0.5, lab, transform=a.transAxes, fontsize=12,
               fontweight="bold", va="center", ha="left",
               bbox=dict(boxstyle="round,pad=0.15", fc="w", ec="none", alpha=0.7))

    fig.tight_layout()
    _save("fig_2_rmse_bias_qq")
    plt.close()

if __name__ == "__main__":
    torch.set_float32_matmul_precision('high')
    t_start = time.perf_counter()
    run_section_9_3()
    print(f"Total time {time.perf_counter() - t_start:.1f} s")
