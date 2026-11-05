import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.raser2 import linear_threshold, offsets, simulate_raser
from src.style import PLUM, use_girlypop

use_girlypop()
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

T1, T2, tau_rd, fwhm, n = 5.0, 3.0, 0.2, 0.3, 101
T2_star = 1 / (1/ T2 + np.pi * fwhm)
simple_rule = -tau_rd / T2_star
w = np.ones(n) / n
print(f"T2* from the linewidth: {T2_star:.2f} s, simple-rule threshold {simple_rule:.3f}")

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))
for shape in ("lorentzian", "gaussian"):
    d = offsets(shape, fwhm, n)
    th = linear_threshold(d, T2, tau_rd)
    print(f"{shape:10s}: exact threshold {th:.3f}")

    grid = np.linspace(-1, 1, 400)
    if shape == "lorentzian":
        dens = (fwhm / 2) / np.pi / (grid**2 + (fwhm / 2) ** 2)
    else:
        sig = fwhm / (2 * np.sqrt(2 * np.log(2)))
        dens = np.exp(-grid**2 / (2 * sig**2)) / (sig * np.sqrt(2 * np.pi))
    ax1.plot(grid, dens, label=shape)

    t_fid = np.linspace(0, 6, 3000)
    fid = np.abs((w[:, None] * np.exp(1j * d[:, None] * t_fid)).sum(axis=0)) * np.exp(-t_fid / T2)
    ax2.plot(t_fid, fid, label=shape)

    pumps = np.linspace(-0.05, -0.6, 12)
    mean_amp, spread = [], []
    for p in pumps:
        t, Mx, My, Mz = simulate_raser(d, 300.0, T1, T2, tau_rd, p * w, n_points=3000,
                                       tilt=1e-5, Mz_start=p * w)
        amp = np.hypot(Mx.sum(0), My.sum(0))[t > 250]
        mean_amp.append(amp.mean())
        spread.append(amp.std())
    line = ax3.errorbar(-pumps, mean_amp, yerr=spread, fmt="o-", ms=4, capsize=2, label=shape)
    ax3.axvline(-th, color=line[0].get_color(), ls=":", lw=1.5)


ax1.set(xlabel="Frequency offset (Hz)", ylabel="Spins per Hz", title="Two lines, same width (0.3 Hz)")
ax1.legend(fontsize=8)
ax2.plot(t_fid, np.exp(-t_fid / T2_star), "--", color=PLUM, lw=1.2, label="exp(−t/T2*)")
ax2.set(xlabel="Time (s)", ylabel="|M transverse|", title="FID: only the Lorentzian is exponential")
ax2.legend(fontsize=8)
ax3.axvline(-simple_rule, color=PLUM, ls="--", lw=1.2, label="simple rule −τ_rd/T2*")
ax3.plot([], [], ":", color=PLUM, label="exact threshold (dotted, per shape)")
ax3.set(xlabel="Pumping strength (−Mz_eq)", ylabel="Mean signal (± fluctuation)",
        title="RASER threshold: line shape matters")
ax3.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "06_inhomogeneity.png", dpi=150)