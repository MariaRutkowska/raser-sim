import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.bloch import simulate
from src.style import PLUM, use_girlypop

use_girlypop()
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

T1, T2, tau_rd = 5.0, 1.0, 0.2
threshold = -tau_rd / T2          # pumping needed to start oscillating
print(f"Threshold: Mz_eq must be below {threshold:.2f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))


t, (Mx, My, Mz) = simulate([1e-3, 0.0, -1.0], 60.0, 0.0, T1, T2,
                           n_points=6000, tau_rd=tau_rd, Mz_eq=-1.0)
ax1.plot(t, np.hypot(Mx, My), label="|M transverse| (the signal)")
ax1.plot(t, Mz, label="Mz")
ax1.axhline(threshold, color=PLUM, ls="--", lw=1, label=f"threshold Mz = {threshold:.1f}")
ax1.set(xlabel="Time (s)", ylabel="Magnetization", title="RASER switches on by itself")
ax1.legend(fontsize=8, loc="lower right")
print(f"Steady signal amplitude: {np.hypot(Mx[-1], My[-1]):.3f}")


pumps = np.linspace(0.0, -1.0, 41)
amps = []
for p in pumps:
    t, (Mx, My, Mz) = simulate([1e-3, 0.0, p], 300.0, 0.0, T1, T2,
                               n_points=3000, tau_rd=tau_rd, Mz_eq=p)
    amps.append(np.hypot(Mx[-1], My[-1]))
theory = np.sqrt(np.clip(tau_rd * (threshold - pumps) / T1, 0, None))
ax2.plot(-pumps, amps, "o", ms=4, label="simulation")
ax2.plot(-pumps, theory, color=PLUM, lw=1.2, label="theory")
ax2.axvline(-threshold, color=PLUM, ls=":", lw=1)
ax2.set(xlabel="Pumping strength  (−Mz_eq)", ylabel="Steady signal amplitude",
        title="Like a laser: nothing below threshold")
ax2.legend(fontsize=8)

fig.tight_layout()
fig.savefig(FIG / "03_raser.png", dpi=150)