import glob
import os
import subprocess
import time
from Bio import SeqIO

in_dir = "../data/gisaid/"
out_dir = "../processed_data/phylogenetic/msa/"
os.makedirs(out_dir, exist_ok=True)


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
	cmd = f"iqtree -s {fin} -pre {out_prefix}"
	output = subprocess.run(cmd, shell=True, capture_output=True)
	duration = time.time() - start 
	print(f"Duration: {str(duration)}")

make_input(in_dir, out_dir)
msa(f"{out_dir}/gisaid.fasta", f"{out_dir}/gisaid.ali")
construct_tree(f"{out_dir}/gisaid.ali", f"{out_dir}/gisaid")