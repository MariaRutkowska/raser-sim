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

T1=T2= 3.37
tau_rd = 0.020
delta = 0.0

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(11,4))

for label, tau in [("no radiation damping", np.inf), (f"tau_rd = {tau_rd * 1e3:.0f} ms", tau_rd)]:
    t, (Mx, My, Mz) = simulate([1.0, 0.0, 0.0], 0.3, delta, T1, T2, tau_rd=tau)
    ax1.plot(t * 1e3, np.hypot(Mx, My), label=f"|M transverse|, {label}")
ax1.plot(t*1e3, Mz, "--", color=PLUM, lw=1.2, label="Mz, with radiation damping")
ax1.set(xlabel="Time (ms)", ylabel="Magnetization", title="After a 90 pulse: signal dies fast")
ax1.legend(fontsize=8, loc="center right")

tilt = np.deg2rad(0.1)
t, (Mx, My, Mz) = simulate([np.sin(tilt), 0.0, -np.cos(tilt)], 0.4, delta, T1, T2, tau_rd=tau_rd)
ax2.plot(t * 1e3, np.hypot(Mx, My), label="|M transverse|")
ax2.plot(t * 1e3, Mz, label="Mz")
ax2.axhline(0, color=PLUM, lw=0.8)
ax2.set(xlabel="Time (ms)", ylabel="Magnetization", title="After a 180° pulse: spins flip back by themselves")
ax2.legend(fontsize=8)


fig.tight_layout()
fig.savefig(FIG / "02_radiation_damping.png", dpi=150)

t_peak = t[np.argmax(np.hypot(Mx, My))]
print(f"Inverted spins: signal peaks after {t_peak * 1e3:.0f} ms")
print(f"Prediction tau_rd * ln(2/tilt) = {tau_rd * np.log(2 / tilt) * 1e3:.0f} ms")

