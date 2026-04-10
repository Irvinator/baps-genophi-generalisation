#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


# Cross-validated out-of-fold confusion matrix
cm = np.array([
    [321, 41],
    [101, 11]
], dtype=int)

# Row-wise percentages: within each true class
row_sums = cm.sum(axis=1, keepdims=True)
pct = cm / row_sums * 100

labels = np.array([
    [f"{cm[0,0]}\n({pct[0,0]:.1f}%)", f"{cm[0,1]}\n({pct[0,1]:.1f}%)"],
    [f"{cm[1,0]}\n({pct[1,0]:.1f}%)", f"{cm[1,1]}\n({pct[1,1]:.1f}%)"],
])

fig, ax = plt.subplots(figsize=(6.6, 5.6), dpi=300)

im = ax.imshow(cm, interpolation="nearest", cmap="Blues")

# Axis ticks and labels
ax.set_xticks([0, 1])
ax.set_yticks([0, 1])
ax.set_xticklabels(["Predicted\nNon-interaction", "Predicted\nInteraction"], fontsize=11)
ax.set_yticklabels(["Actual\nNon-interaction", "Actual\nInteraction"], fontsize=11)

# Cell annotations
threshold = cm.max() / 2.0
for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        ax.text(
            j, i, labels[i, j],
            ha="center", va="center",
            color="white" if cm[i, j] > threshold else "black",
            fontsize=12
        )

ax.set_title("Cross-validated confusion matrix", fontsize=13, pad=10)
ax.set_xlabel("Predicted class", fontsize=11)
ax.set_ylabel("True class", fontsize=11)

# Clean up layout
for spine in ax.spines.values():
    spine.set_visible(True)

plt.tight_layout()

out = "/mnt/d/baps-genophi-generalisation/report_figures/cv_confusion_matrix_final.png"
plt.savefig(out, dpi=300, bbox_inches="tight")
plt.savefig("/mnt/d/baps-genophi-generalisation/report_figures/cv_confusion_matrix_final.pdf", bbox_inches="tight")
print(f"Saved to: {out}")
plt.show()
