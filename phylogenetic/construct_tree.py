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


def msa(mode, fasta_fn, align_fn, existing_alignment=None):
	assert mode in ["denovo", "add"], \
		"Argument mode must be denovo or add."

	print("Aligning with MAFFT.")
	start = time.time()

	if mode == "denovo":
		print("mode = denovo.")
		cmd = f"mafft {fasta_fn} > {align_fn}"
	else:
		assert existing_alignment is not None, \
			"Argument existing_alignment must be provided."
		print("mode = add.")
		cmd = f"mafft --add {fasta_fn} {existing_alignment} > {align_fn}"

	print(f"Command: {cmd}")
	output = subprocess.run(cmd, shell=True, capture_output=True)
	
	duration = time.time() - start
	print("Finished MSA.")
	print(f"Time lapsed: {str(duration)}")


def construct_tree(software, msa_fn, out_prefix=None,\
	out_fn=None, log_fn=None):
	assert software in ["iqtree", "FastTree"], \
		"Software must be iqtree or FastTree."

	print("Contructing evolutionary tree.")
	start = time.time()

	if software == "iqtree":
		print("Software = iqtree.")
		assert out_prefix is not None, \
			"iqtree requires the out_prefix argument." 
		cmd = f"iqtree -s {msa_fn} -pre {out_prefix} -m GTR"
		if os.path.exists(f"{out_prefix}.ckp.gz"):
			cmd += " -redo"
	elif software == "FastTree":
		print("Software = FastTree.")
		assert (out_fn is not None) and (log_fn is not None), \
			"FastTree requires out_fn and log_fn arguments."
		cmd = f"FastTree -nt -gtr -out {out_fn} -log {log_fn} {msa_fn}"
	else:
		pass

	print(f"Command: {cmd}")
	output = subprocess.run(cmd, shell=True, capture_output=True)

	duration = time.time() - start 
	print("Finished tree construction.")
	print(f"Time lapsed: {str(duration)}")



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