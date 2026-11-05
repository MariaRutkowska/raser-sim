import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.raser2 import linear_threshold, offsets
from src.style import PLUM, use_girlypop

use_girlypop()
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

T2, tau_rd, fwhm = 3.0, 0.2, 0.3
T2_star = 1 / (1 / T2 + np.pi * fwhm)
ns = [5, 11, 21, 51, 101, 201, 501, 1001, 10001]

print("| Groups | Threshold, Lorentzian | Threshold, Gaussian |")
print("|---|---|---|")
results = {"lorentzian": [], "gaussian": []}
for n in ns:
    row = []
    for shape in results:
        th = linear_threshold(offsets(shape, fwhm, n), T2, tau_rd)
        results[shape].append(th)
        row.append(f"{th:.4f}")
    print(f"| {n} | {row[0]} | {row[1]} |")

fig, ax = plt.subplots(figsize=(7, 4))
for shape, ths in results.items():
    ax.semilogx(ns, -np.array(ths), "o-", ms=5, label=shape)
ax.axhline(tau_rd / T2_star, color=PLUM, ls="--", lw=1.2, label="simple rule τ_rd/T2*")
ax.axvline(101, color=PLUM, ls=":", lw=1, label="101 groups used in the simulations")
ax.set(xlabel="Number of spin groups", ylabel="Threshold pumping (−Mz_eq)",
       title="The result stops changing after ~20 groups")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIG / "07_convergence.png", dpi=150)