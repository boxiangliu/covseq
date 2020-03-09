import subprocess
import os
import shutil
import glob
import click

def run_glimmer3_iterated(fasta_fn, out_dir):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	cmd = f"../src/glimmer3.02/scripts/g3-iterated.csh {fasta_fn} {out_dir}/iterated"
	cmd = cmd.split()
	print(cmd)
	output = subprocess.run(cmd, capture_output=True)
	return output

def clean(in_dir, out_dir):
	fout = None
	with open(f"{in_dir}/iterated.predict", "r") as fin:
		for line in fin:
			if line.startswith(">"):
				strain = line.strip().split()[0].replace(">", "")
				if fout:
					fout.close()
				if not os.path.exists(f"{out_dir}/{strain}"):
					os.makedirs(f"{out_dir}/{strain}")
				fout = open(f"{out_dir}/{strain}/{strain}.predict", "w")
			fout.write(line)
	for i in glob.glob(f"{in_dir}/iterated.*"):
		os.remove(i)
	os.rmdir(in_dir)


@click.command()
@click.option("-f", "--fasta", type=str, required=True, \
	help="Fasta file containing one or more virus strains.")
@click.option("-o", "--out_dir", type=str, required=False, \
	help="Output directory", default="results", show_default=True)
def predict(fasta_fn, out_dir):
	run_glimmer3_iterated(fasta_fn, f"{out_dir}/orf_prediction/")
	clean(f"{out_dir}/orf_prediction/", out_dir)

