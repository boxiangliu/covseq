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


msa(mode="related", \
	fasta_fn="../processed_data/phylogenetic/filter_distant_seq/keep.fasta", \
	align_fn="../processed_data/phylogenetic/construct_tree/keep.ali", \
	ref_fasta_fn="data/NC_045512.2.fasta", uppercase=True)

msa(mode="related", \
	fasta_fn="../processed_data/phylogenetic/append_new_seq/keep_and_new.fasta", \
	align_fn="../processed_data/phylogenetic/construct_tree/keep_and_new.ali", \
	ref_fasta_fn="data/NC_045512.2.fasta", uppercase=True)


def msa(mode, fasta_fn, align_fn, existing_alignment=None, ref_fasta_fn=None, uppercase=False):
	assert mode in ["denovo", "add", "related"], \
		"Argument mode must be denovo or add or related."

	print("Aligning with MAFFT.")
	start = time.time()

	if mode == "denovo":
		print("mode = denovo.")
		cmd = f"mafft --thread -1 {fasta_fn} > {align_fn}"

	elif mode == "add":
		assert existing_alignment is not None, \
			"existing_alignment must be provided."
		print("mode = add.")
		cmd = f"mafft --thread -1 --add {fasta_fn} {existing_alignment} > {align_fn}"

	elif mode == "related":
		assert ref_fasta_fn is not None, \
			"ref_fasta_fn must be provided."
		print("mode = related (input sequence must have % identity ~ 99%.")
		cmd = f"mafft --auto --thread -1 --addfragments {fasta_fn} {ref_fasta_fn} > {align_fn}"
	else:
		pass
	print(f"Command: {cmd}")
	output = subprocess.run(cmd, shell=True, capture_output=True)
	
	msa_time = time.time()
	duration = msa_time - start
	print(f"MSA time: {str(duration)}")

	if uppercase:
		print("Convert to uppercase")
		with open(align_fn + ".up", "w") as f:
			for fasta in SeqIO.parse(align_fn, "fasta"):
				fasta.seq = fasta.seq.upper()
				SeqIO.write(fasta, f, "fasta")
		os.remove(align_fn)
		os.rename(align_fn + ".up", align_fn)


options = "-threads 12 -double-precision -ext AVX2 -fastexp 2 "
output = construct_tree(software="VeryFastTree", \
	msa_fn="../processed_data/phylogenetic/construct_tree/keep.ali", \
	out_fn="../processed_data/phylogenetic/construct_tree/keep.vft.nh", \
	log_fn="../processed_data/phylogenetic/construct_tree/keep.vft.log", \
	options=options)

options = "-threads 12 -double-precision -ext AVX2 -fastexp 2 -fastest "
output = construct_tree(software="VeryFastTree", \
	msa_fn="../processed_data/phylogenetic/construct_tree/keep.ali", \
	out_fn="../processed_data/phylogenetic/construct_tree/keep.fastest.vft.nh", \
	log_fn="../processed_data/phylogenetic/construct_tree/keep.fastest.vft.log", \
	options=options)


options = "-threads 12 -double-precision -ext AVX2 -fastexp 2 "
output = construct_tree(software="VeryFastTree", \
	msa_fn="../processed_data/phylogenetic/construct_tree/keep_and_new.ali", \
	out_fn="../processed_data/phylogenetic/construct_tree/keep_and_new.vft.nh", \
	log_fn="../processed_data/phylogenetic/construct_tree/keep_and_new.vft.log", \
	options=options)

def construct_tree(software, msa_fn, out_prefix=None,\
	out_fn=None, log_fn=None, options=""):
	assert software in ["iqtree", "FastTree", "VeryFastTree"], \
		"Software must be iqtree or FastTree or VeryFastTree."

	print("Contruct evolutionary tree.")
	start = time.time()

	if software == "iqtree":
		print("Software = iqtree.")
		assert out_prefix is not None, \
			"iqtree requires the out_prefix argument." 
		cmd = f"iqtree -s {msa_fn} -pre {out_prefix} -m GTR"
		if os.path.exists(f"{out_prefix}.ckp.gz"):
			cmd += " -redo"

	elif software == "FastTree" or software == "VeryFastTree":
		assert (out_fn is not None) and (log_fn is not None), \
			"FastTree requires out_fn and log_fn arguments."
		print(f"Software = {software}")
		cmd = f"{software} -nt -gtr -out {out_fn} -log {log_fn}"
		if options != "":
			cmd += " " + options
		cmd += " " + msa_fn

	else:
		pass

	print(f"Command: {cmd}")
	output = subprocess.run(cmd, shell=True, capture_output=True)

	duration = time.time() - start 
	print("Finished tree construction.")
	print(f"Time lapsed: {str(duration)}")
	return output



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