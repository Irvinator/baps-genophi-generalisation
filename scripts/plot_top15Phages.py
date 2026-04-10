import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("/mnt/d/baps-genophi-generalisation/aim2_ecoli96/aim2_ecoli96_ranked_candidates.csv")

top = df.head(15)

plt.figure()
plt.barh(top["phage_id"][::-1], top["candidate_priority_score"][::-1])
plt.xlabel("Candidate priority score")
plt.ylabel("Phage ID")
plt.title("Top-ranked candidate phages for TXTL rebootability screening")
plt.tight_layout()

plt.savefig("/mnt/d/baps-genophi-generalisation/aim2_ecoli96/aim2_top15_barplot.png", dpi=300)
plt.show()
