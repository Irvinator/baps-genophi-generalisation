import os
import pandas as pd
from Bio import SeqIO

# input directory
input_dir = "/mnt/d/baps-genophi-generalisation/aim2_ecoli96/millard_phage_genomes"

rows = []

for file in os.listdir(input_dir):
    if file.endswith(".fa"):
        filepath = os.path.join(input_dir, file)
        
        total_length = 0
        gc_count = 0
        contigs = 0
        
        for record in SeqIO.parse(filepath, "fasta"):
            seq = str(record.seq).upper()
            total_length += len(seq)
            gc_count += seq.count("G") + seq.count("C")
            contigs += 1
        
        gc_content = gc_count / total_length if total_length > 0 else 0
        
        rows.append({
            "phage_id": file.replace(".fa", ""),
            "genome_length": total_length,
            "gc_content": gc_content,
            "num_contigs": contigs
        })

df = pd.DataFrame(rows)

output_path = "/mnt/d/baps-genophi-generalisation/aim2_ecoli96/aim2_ecoli96_basic_features.csv"
df.to_csv(output_path, index=False)

print("Saved:", output_path)
print(df.head())
