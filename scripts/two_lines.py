import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.raser2 import simulate_raser
from src.style import PLUM, use_girlypop

use_girlypop()
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

T1, T2, tau_rd, pump = 5.0, 1.0, 0.2, -1.0
regimes = {"far apart (2 Hz): two lines, beating": 2.0,
           "close (0.2 Hz): chaos": 0.2,
           "very close (0.1 Hz): lines lock into one": 0.1}


fig, axes = plt.subplots(2, 3, figsize=(13, 6))
for col, (title, gap) in enumerate(regimes.items()):
    deltas = 2 * np.pi * np.array([gap / 2, -gap / 2])
    t, Mx, My, Mz = simulate_raser(deltas, 400.0, T1, T2, tau_rd, [pump, pump], n_points=40000)
    signal = Mx.sum(axis=0) + 1j * My.sum(axis=0)   

    late = t > 200                               
    ax = axes[0, col]
    ax.plot(t[late], np.abs(signal[late]), lw=1)
    ax.set(xlim=(200, 240), ylim=(0, 1.2 * np.abs(signal[late]).max()),
           xlabel="Time (s)", ylabel="|signal|", title=title)

    s, tt = signal[late], t[late]
    spec = np.abs(np.fft.fftshift(np.fft.fft(s * np.hanning(len(s)))))
    freq = np.fft.fftshift(np.fft.fftfreq(len(s), tt[1] - tt[0]))
    ax = axes[1, col]
    ax.semilogy(freq, spec / spec.max(), lw=1)
    for d in (gap / 2, -gap / 2):
        ax.axvline(d, color=PLUM, ls=":", lw=1)  
    ax.set(xlim=(-2.5, 2.5), ylim=(1e-5, 2), xlabel="Frequency (Hz)", ylabel="Spectrum")
fig.tight_layout()
fig.savefig(FIG / "04_two_lines.png", dpi=150)


fig, ax = plt.subplots(figsize=(7, 4))
for title, gap in [("2 Hz (regular)", 2.0), ("0.2 Hz (chaos)", 0.2)]:
    deltas = 2 * np.pi * np.array([gap / 2, -gap / 2])
    runs = []
    for tilt in (1e-3, 1e-3 + 1e-9):
        t, Mx, My, Mz = simulate_raser(deltas, 500.0, T1, T2, tau_rd, [pump, pump],
                                       n_points=50000, tilt=tilt)
        runs.append(np.hypot(Mx.sum(axis=0), My.sum(axis=0)))
    ax.semilogy(t, np.abs(runs[0] - runs[1]) + 1e-16, lw=0.8, label=title)
ax.set(xlabel="Time (s)", ylabel="Difference between two runs",
       title="Chaos: a 1e-9 change in the start grows to everything")
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(FIG / "04_chaos.png", dpi=150)
print("Figures saved")