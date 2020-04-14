import glob
import os
import subprocess
import time
from Bio import SeqIO
import click

def make_input(in_dir, out_dir):
	with open(f"{out_dir}/gisaid.fasta", "w") as fout:
		counter = 0
		for fn in glob.glob(f"{in_dir}/*.fasta"):
			fasta = next(SeqIO.parse(fn, "fasta"))
			fasta.id = fasta.description.split("|")[1]
			if len(fasta.seq) < 29000:
				print(f"{fasta.id} only has {len(fasta.seq)} nucleotides. Skipping...")
			else:
				counter += 1
				fout.write(f">{fasta.id}\n")
				fout.write(str(fasta.seq) + "\n")
		print(f"Total: {counter} sequences.")


def msa(fin, fout):
	print("Aligning with MAFFT...")
	start = time.time()
	cmd = f"mafft {fin} > {fout}"
	output = subprocess.run(cmd, shell=True, capture_output=True)
	duration = time.time() - start
	print(f"Duration: {str(duration)}")


def construct_tree(fin, out_prefix):
	print("Contructing evolutionary tree...")
	start = time.time()
	cmd = f"iqtree -s {fin} -pre {out_prefix} -m GTR -T AUTO"
	output = subprocess.run(cmd, shell=True, capture_output=True)
	duration = time.time() - start 
	print(f"Duration: {str(duration)}")


@click.command()
@click.option("--in_fn", "-i", type=str, help="Input FASTA file.")
@click.option("--out_dir", "-o", type=str, help="Output directory.")
@click.option("--make_tree", "-t", is_flag=True, default=False)
def main(in_fn, out_dir, make_tree):
	print("##############")
	print("# MSA / Tree #")
	print("##############")
	print(f"FASTA: {in_fn}")
	print(f"Output: {out_dir}")

	os.makedirs(out_dir, exist_ok=True)
	msa(in_fn, f"{out_dir}/preprocessed.ali")
	if make_tree: construct_tree(f"{out_dir}/preprocessed.ali", f"{out_dir}/preprocessed")


if __name__ == "__main__":
	main()