import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("/home/irvin/data/BAPS/outputs/phage_cluster_overlap.csv")

plt.hist(df["mapped_proteins"], bins=40)

plt.xlabel("Proteins mapped to GenoPHI clusters")
plt.ylabel("Number of BAPS phages")
plt.title("Protein family overlap between BAPS phages and GenoPHI training phages")

plt.savefig("phage_cluster_overlap.png", dpi=300)
