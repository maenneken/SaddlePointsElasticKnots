"""
Cumulative success rate plot for motion planning benchmarks.

Usage:
    Fill in the `datasets` dictionary below with your own iteration counts
    for each knot pair / configuration you want to compare, then run the
    script. It produces both a .pdf (for LaTeX) and a .png (for quick
    viewing) in the same folder.
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

# ---- Style settings (LaTeX-like serif font, consistent with thesis) ----
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 11
matplotlib.rcParams['axes.labelsize'] = 12
matplotlib.rcParams['axes.titlesize'] = 12
matplotlib.rcParams['legend.fontsize'] = 10

# ---- 1. Insert your data here ----
# Each entry: "label": [list of iteration counts for successful runs]
datasets = {
    r"$4_1$ (easy)": [59, 69, 90, 114, 114, 130, 144, 152, 153, 163, 164,
                       173, 188, 189, 197, 199, 204, 205, 211, 212, 216,
                       217, 225, 233, 247, 247, 255, 271, 285, 285, 291,
                       297, 303, 305, 319, 345, 346, 423, 435, 449, 457,
                       476, 488, 501, 514, 517, 536, 544, 546, 666],
    r"$4_1$ (medium)": [378, 413, 450, 666, 677, 690, 792, 904, 1257, 1330, 1553, 1600, 1824, 1834, 1871, 2220, 2579, 2655, 2880, 3303, 3532, 3949, 4516, 4707, 5095, 5337, 5696, 5888, 6084, 6310, 7257, 7931, 8614, 10771, 11158, 11486, 12937, 13057, 15310, 15428, 16989, 18104, 19259, 20047, 20166, 21079, 21169, 25568, 36066, 43649],
    r"$5_2$ (hard)":   [2331, 3698, 4017, 4480, 6483, 6625, 7731, 7962, 8134, 8774, 8940, 9599, 10234, 10280, 10424, 11194, 11240, 11380, 11799, 13990, 14078, 14224, 14687, 16020, 16716, 17342, 17948, 18009, 18977, 19736, 20705, 22588, 24795, 25443, 25681, 26952, 28086, 28162, 28311, 30278, 30976, 34172, 41582, 46176, 46392,  np.inf , np.inf, np.inf, np.inf, np.inf]
}

# ---- 2. Colors (extend if you have more than 3 datasets) ----
colors = ['#1f5e8c', '#c0392b', '#27964e', '#8e44ad']

# ---- 3. Plotting ----
fig, ax = plt.subplots(figsize=(6.5, 4.2))

for (label, data), color in zip(datasets.items(), colors):
    data = sorted(data)
    n = len(data)
    x = [0] + data
    y = [0] + list(np.arange(1, n + 1) / n * 100)
    ax.step(x, y, where='post', color=color, linewidth=2, label=label)
    ax.scatter(data, np.arange(1, n + 1) / n * 100, color=color, s=14, zorder=3)

ax.set_xlabel('Number of iterations')
ax.set_ylabel('Success rate (\\%)')
ax.set_ylim(0, 105)
ax.set_xlim(left=0)
ax.grid(True, linestyle='--', alpha=0.4)
ax.legend(loc='lower right', frameon=False)

fig.tight_layout()
fig.savefig('success_rate_curve.pdf', bbox_inches='tight')
fig.savefig('success_rate_curve.png', dpi=200, bbox_inches='tight')
print("Saved success_rate_curve.pdf and .png")
