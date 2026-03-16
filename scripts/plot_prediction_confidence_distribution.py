import pandas as pd
import matplotlib.pyplot as plt

# Load prediction results
df = pd.read_csv("/home/irvin/projects/GenoPHI/results/baps_predict/strain_median_predictions.csv")

# Plot confidence distribution
plt.figure(figsize=(8,6))

plt.hist(df["Confidence"], bins=40)

plt.xlabel("Predicted interaction probability")
plt.ylabel("Number of strain–phage pairs")
plt.title("Distribution of GenoPHI predicted interaction probabilities for BAPS dataset")

plt.tight_layout()

plt.savefig("prediction_confidence_distribution.png", dpi=300)

print("Saved prediction_confidence_distribution.png")
