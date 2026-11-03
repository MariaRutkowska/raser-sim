"""IMPORTANT: free precession of protons in pure water, and the spectrum.

T1 of degassed water at 25 C: 3.37 s (Simpson & Carr, Phys. Rev. 111, 1201, 1958).
Small, fast-tumbling molecule (extreme narrowing): T2 = T1.
In a real magnet the field is not perfectly uniform, so the signal decays
faster, with T2* instead of T2. Here T2* = 0.3 s as an example.
"""

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

delta = 2 * np.pi * 5.0
T1 = 3.37
cases = { "ideal magnet: T2 = T1 = 3.37 s": T1,
         "real magnet: T2* = 0.3 s": 0.3}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11,4))
for label, T2 in cases.items():
    t, (Mx, My, Mz) = simulate([1.0, 0.0, 0.0], t_end=20.0, delta=delta, T1=T1, T2=T2, n_points=8000)
    ax1.plot(t, Mx, lw=1.2, label=label)

    signal = Mx +1j * My
    dt = t[1] - t[0]
    spec = np.fft.fftshift(np.fft.fft(signal))
    freq = np.fft.fftshift(np.fft.fftfreq(len(t), dt))
    ax2.plot(freq, np.abs(spec) / np.abs(spec).max(), label=label)
    print(f"{label}: linewidth 1/(pi*T2) = {1 / (np.pi * T2):.3f} Hz")

ax1.plot(t, Mz, "--", color=PLUM, lw=1.2, label="Mz (recovers with T1)")
ax1.set(xlim=(0, 10), xlabel="Time (s)", ylabel="Magnetization", title="Water protons: FID")
ax1.legend(fontsize=8, loc="lower right")
ax2.set(xlim=(2, 8), xlabel="Frequency (Hz)", ylabel="Normalised |spectrum|", title="Spectrum: shorter T2 = wider line")
fig.tight_layout()
fig.savefig(FIG / "01_free_precession.png", dpi=150)