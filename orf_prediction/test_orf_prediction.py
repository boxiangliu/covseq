import pytest
from orf_prediction import *

def test_run_glimmer3_iterated():
	fasta_fn = "example.fasta"
	out_dir = "results/example/"
	print("a", flush=True)
	output = run_glimmer3_iterated(fasta_fn, out_dir)
	assert "Command not found" not in output.stderr.decode("utf-8")


def test_clean():
	in_dir = "results/example/"
	out_dir = "results/annotation/"
	clean(in_dir, out_dir)
	assert os.path.exists(f"{out_dir}/NC_045512.2/NC_045512.2.predict")
	assert not os.path.exists(in_dir)


def test_predict():
	fasta_fn = "example.fasta"
	out_dir = "results/annotation/"
	predict(fasta_fn, out_dir)
	assert os.path.exists(f"{out_dir}/NC_045512.2/NC_045512.2.predict")
	assert os.path.exists(f"{out_dir}/EPI_ISL_408669/EPI_ISL_408669.predict")
