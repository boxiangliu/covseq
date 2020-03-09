import pytest
from orf_prediction import *

def test_run_glimmer3_iterated():
	fasta_fn = "example.fasta"
	out_dir = "results/example/"
	output = run_glimmer3_iterated(fasta_fn, out_dir)
	assert output.returncode == 0