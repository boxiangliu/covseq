import glob
import gzip
from Bio import SeqIO
in_dir = "../data/"
sub_dirs = ["cngb", "embl", "gisaid", "ncbi"]
out_fn = "../data/aggregated/fasta/raw.fasta"

print("#########################")
print("# Concatenate all FASTA #")
print("#########################")
print(f"Input: {in_dir}")
print(f"Output: {out_fn}")
counter = 0
with open(out_fn, "w") as fout:
	for sd in sub_dirs:
		wd = f"{in_dir}/{sd}/fasta/"
		print(f"Working on directory {wd}")
		for fn in glob.glob(f"{wd}/*.fasta"):
			for record in SeqIO.parse(fn, "fasta"):
				SeqIO.write(record, fout, "fasta")
				counter += 1
print(f"Combined {counter} sequences!")

