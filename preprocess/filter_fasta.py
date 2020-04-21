import os
import hashlib
import click
from Bio import SeqIO


def remove_duplicates(in_fn, out_dir):
	print("Removing duplicates.")
	hash_set = set()
	out_table_fn = f"{out_dir}/duplicated.tsv"
	out_fasta_fn = f"{out_dir}/non_redundant.fasta"
	with open(out_table_fn, "w") as out_table, \
		open(out_fasta_fn, "w") as out_fasta:
		out_table.write("ID\tduplicated\n")

		for record in SeqIO.parse(in_fn, "fasta"):
			sha256 = hashlib.sha256(str(record.seq).encode("utf-8")).hexdigest()
			if sha256 not in hash_set:
				hash_set.add(sha256)
				out_table.write(f"{record.description}\t0\n")
				SeqIO.write(record, out_fasta, "fasta")
			else:
				out_table.write(f"{record.description}\t1\n")
	return out_table_fn, out_fasta_fn


def remove_duplicates_by_ID(in_fn, out_dir):
	print("Removing duplicates.")
	hash_set = set()
	out_table_fn = f"{out_dir}/duplicated.tsv"
	out_fasta_fn = f"{out_dir}/non_redundant.fasta"
	with open(out_table_fn, "w") as out_table, \
		open(out_fasta_fn, "w") as out_fasta:
		out_table.write("ID\tduplicated\n")

		for record in SeqIO.parse(in_fn, "fasta"):
			if record.id not in hash_set:
				hash_set.add(record.id)
				out_table.write(f"{record.id}\tNO\n")
				SeqIO.write(record, out_fasta, "fasta")
			else:
				out_table.write(f"{record.id}\tYES\n")
	return out_table_fn, out_fasta_fn



def get_genome_length(in_fn, out_fn):
	print("Getting sequence lengths.")
	with open(out_fn, "w") as fout:
		for record in SeqIO.parse(in_fn, "fasta"):
			length = len(record.seq)
			fout.write(f"{record.description}\t{length}\n")


def filter_complete_genome(in_fn, out_fn, cutoff=25000):
	print("Filtering for complete genomes.")
	with open(out_fn, "w") as fout:
		for record in SeqIO.parse(in_fn, "fasta"):
			length = len(record.seq)
			if length >= cutoff:
				SeqIO.write(record, fout, "fasta")


def count_ambiguous_base(in_fn, out_fn):
	print("Counting ambiguous bases.")
	with open(out_fn, "w") as fout:
		for record in SeqIO.parse(in_fn, "fasta"):
			ambi_base = 0
			length = len(record.seq)
			for b in str(record.seq):
				if b.upper() not in "ATCGU":
					ambi_base += 1
			proportion = ambi_base / length
			fout.write(f"{record.description}\t{ambi_base}\t{length}\t{proportion}\n")


@click.command()
@click.option("--in_fn", "-i", type=str, help="Input file.")
@click.option("--out_dir", type=str, help="Directory to put intermediate files.")
@click.option("--final_fn", type=str, help="Final fasta file.")
def main(in_fn, out_dir, final_fn):
	print("#############################")
	print("# Filtering FASTA sequences #")
	print("#############################")
	print(f"Input: {in_fn}")
	print(f"Output: {final_fn}")
	
	os.makedirs(out_dir, exist_ok=True)
	_, non_redundant = remove_duplicates_by_ID(in_fn, out_dir)
	get_genome_length(non_redundant, f"{out_dir}/genome_lengths.tsv")
	filter_complete_genome(non_redundant, final_fn, cutoff=25000)
	count_ambiguous_base(final_fn, f"{out_dir}/ambiguous_bases.tsv")


if __name__ == "__main__":
	main()