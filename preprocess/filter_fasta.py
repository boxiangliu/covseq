import os
from Bio import SeqIO

in_fn = "../data/aggregated/fasta/raw.fasta"
out_dir = "../processed_data/preprocess/filter_fasta/"
os.makedirs(out_dir, exist_ok=True)
out_fn = "../data/aggregated/fasta/processed.fasta"

print("#############################")
print("# Filtering FASTA sequences #")
print("#############################")
print(f"Input: {in_fn}")
print(f"Output: {out_fn}")

print("Getting sequence lengths.")
with open(f"{out_dir}/genome_lengths.tsv", "w") as fout:
	for record in SeqIO.parse(in_fn, "fasta"):
		length = len(record.seq)
		fout.write(f"{record.description}\t{length}\n")

print("Filtering for complete genomes.")
cutoff = 25000
with open(f"{out_dir}/complete_genomess.fasta", "w") as fout:
	for record in SeqIO.parse(in_fn, "fasta"):
		length = len(record.seq)
		if length >= cutoff:
			SeqIO.write(record, fout, "fasta")

