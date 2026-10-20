"""Regime diagram of the terminal state in the (beta, theta_s) plane.

Figure 4 of the submitted manuscript. Reference parameters of Table 3, theta_v = 100.
The three regions are the stability domains of P212, of the pinned
equilibrium P51, and of P112, whose boundaries are theta_s = I212(beta)
and theta_s = I112(beta).

The diagram is computed from the analytic stability classification of
Theorem 2 and Proposition 4; the six open circles are the points verified by
direct integration with src/robustness.py
(results/fig7_checkpoints_corrected.csv).
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"pdf.fonttype":42,"ps.fonttype":42})

LAM, PV, PB, GAM, D = 10.0, 0.897, 0.7934, 0.143, 0.007

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "Figure_4.pdf")
os.makedirs(os.path.dirname(OUT), exist_ok=True)


def I112(b):
    a1 = (GAM + D) * (1 - PB) * b ** 2
    b1 = b * (D * (GAM + D) * (2 - PB) - b * LAM * (1 - PB))
    c1 = D * (D * (GAM + D) - b * LAM * (1 - PB * PV))
    return (-b1 + np.sqrt(b1 * b1 - 4 * a1 * c1)) / (2 * a1)


def I212(b):
    return PV * LAM / (GAM + D) - D / ((1 - PB) * b)


betas = np.linspace(0.001, 0.006, 1200)
fig, ax = plt.subplots(figsize=(4.5, 3.3))
# Analytic bands remain vector paths in the submission PDF.
ax.fill_between(betas, 0, I212(betas), color="#cfe3f7", linewidth=0)
ax.fill_between(betas, I212(betas), I112(betas), color="#ffe2b8", linewidth=0)
ax.fill_between(betas, I112(betas), 100, color="#d8ecc8", linewidth=0)
ax.plot(betas, I212(betas), "-", color="#1f3b57", lw=1.3)
ax.plot(betas, I112(betas), "-", color="#1f3b57", lw=1.3)
ax.axvline(0.003, ls=":", lw=1.0, color="#b3261e")
ax.plot([0.003], [I212(0.003)], "o", ms=4, mfc="white", mec="#1f3b57", mew=1.0)
ax.plot([0.003], [I112(0.003)], "o", ms=4, mfc="white", mec="#1f3b57", mew=1.0)

# points checked by direct integration (src/robustness.py)
for bx in (0.002, 0.005):
    for ty in (20.0, 45.0, 70.0) if bx == 0.002 else (30.0, 54.0, 80.0):
        ax.plot([bx], [ty], "o", ms=3.6, mfc="none", mec="#b3261e", mew=0.9)

ax.text(0.00128, 12, r"$P_{212}$", fontsize=9, color="#1f3b57")
ax.text(0.00355, 82, r"$P_{112}$", fontsize=9, color="#1f3b57")
ax.text(0.00415, 56.6, r"pinned, $I^*=\theta_s$", fontsize=7.5, color="#7a4a12",
        rotation=11, ha="center", va="center")
ax.text(0.00515, 50.2, r"$\theta_s=I_{212}$", fontsize=7, color="#1f3b57", ha="center")
ax.text(0.00515, 63.4, r"$\theta_s=I_{112}$", fontsize=7, color="#1f3b57", ha="center")
ax.text(0.00303, 2, r"$\beta=0.003$", fontsize=6.5, rotation=90, color="#b3261e")

ax.set_xlabel(r"$\beta$")
ax.set_ylabel(r"$\theta_s$")
ax.set_xlim(0.001, 0.006)
ax.set_ylim(0, 100)
ax.tick_params(labelsize=7)
plt.tight_layout()
plt.savefig(OUT)
plt.close()
print("saved", OUT)
