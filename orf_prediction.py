import subprocess

def predict(fasta_fn, out_dir):
	run_glimmer3(fasta_fn, out_dir)
	clean(in_dir, out_dir)

def run_glimmer3_iterated(fasta_fn, out_dir):
	cmd = f"../src/glimmer3.02/scripts/g3-iterated.csh {fasta_fn} {out_dir}/iterated"
	cmd = cmd.split()
	output = subprocess.run(cmd, capture_output=True)
	return output

def clean(in_dir, out_dir):
	pass