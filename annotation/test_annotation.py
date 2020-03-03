import pytest
from annotate import *

def test_annotate():
	ids = get_fasta_ids("example.fasta")
	assert ids == ['NC_045512.2', 'EPI_ISL_408669']


def test_create_metadata():
	create_metadata(["A","B"], out_dir="results/example/")
	with open("results/example/metadata.csv", "r") as f:
		line = f.readline()
		line = f.readline()
		assert line.strip().split(",")[0] == "A"
		assert line.strip().split(",")[-1] == "0"
	os.rmdir("results/example/")

def test_call_vapid():
	fasta_fn = "example.fasta"
	metadata_fn = "example_metadata.csv"
	output = call_vapid(fasta_fn, metadata_fn, \
		out_dir="results/example/vapid/")
	assert output.returncode == 0
