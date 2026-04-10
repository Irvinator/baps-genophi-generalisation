import matplotlib.pyplot as plt
import numpy as np

#Prior array 
'''cm = np.array([
    [23, 0],   # Actual Negative
    [121, 0]   # Actual Positive
])'''

#New confusion matrix from CV-5 on final dataset
cm = np.array([[321, 41],
               [101, 11]])

fig, ax = plt.subplots(figsize=(5.5, 4.5))
im = ax.imshow(cm, cmap="Blues")

ax.set_xticks([0, 1])
ax.set_yticks([0, 1])
ax.set_xticklabels(["Predicted Negative", "Predicted Positive"])
ax.set_yticklabels(["Actual Negative", "Actual Positive"])

for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        ax.text(j, i, str(cm[i, j]), ha="center", va="center", color="black", fontsize=12)

ax.set_title("Cross-Validated Confusion Matrix")
plt.tight_layout()
plt.savefig("cv_confusion_matrix_heatmap.png", dpi=300)
print("Saved confusion_matrix_heatmap.png")
