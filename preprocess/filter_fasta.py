import os
import hashlib
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


def remove_duplicates(in_fn, out_dir):
	print("Removing duplicates.")
	hash_set = set()
	with open(f"{out_dir}/duplicated.tsv", "w") as out_table, \
		open(f"{out_dir}/non_redundant.fasta", "w") as out_fasta:
		out_table.write("ID\tduplicated\n")

		for record in SeqIO.parse(in_fn, "fasta"):
			sha256 = hashlib.sha256(str(record.seq).encode("utf-8")).hexdigest()
			if sha256 not in hash_set:
				hash_set.add(sha256)
				out_table.write(f"{record.description}\t0\n")
				SeqIO.write(record, out_fasta, "fasta")
			else:
				out_table.write(f"{record.description}\t1\n")


def get_genome_length(in_fn, out_dir):
	print("Getting sequence lengths.")
	with open(f"{out_dir}/genome_lengths.tsv", "w") as fout:
		for record in SeqIO.parse(in_fn, "fasta"):
			length = len(record.seq)
			fout.write(f"{record.description}\t{length}\n")


def filter_complete_genome(in_fn, out_dir):
	print("Filtering for complete genomes.")
	cutoff = 25000
	with open(f"{out_dir}/complete_genome.fasta", "w") as fout:
		for record in SeqIO.parse(in_fn, "fasta"):
			length = len(record.seq)
			if length >= cutoff:
				SeqIO.write(record, fout, "fasta")


def count_ambiguous_base(in_fn, out_dir):
	print("Count ambiguous bases.")
	with open(f"{out_dir}/ambiguous_bases.tsv", "w") as fout:
		for record in SeqIO.parse(in_fn, "fasta"):
			ambi_base = 0
			length = len(record.seq)
			for b in str(record.seq):
				if b.upper() not in "ATCGU":
					ambi_base += 1
			proportion = ambi_base / length
			fout.write(f"{record.description}\t{ambi_base}\t{length}\t{proportion}\n")


def main():
	remove_duplicates(in_fn, out_dir)
	# get_genome_length(in_fn, out_dir)
	# filter_complete_genome(in_fn, out_dir)
	# count_ambiguous_base(f"{out_dir}/complete_genome.fasta", out_dir)

if __name__ == "__main__":
	main()