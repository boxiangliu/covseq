import pytest
import shutil
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
	shutil.rmtree("results/example/")

def test_call_vapid():
	fasta_fn = "example.fasta"
	metadata_fn = "example_metadata.csv"
	output = call_vapid(fasta_fn, metadata_fn, \
		out_dir="results/example/vapid/")
	assert output.returncode == 0

def test_parse_vapid():
	out_dir = "results/example/vapid/"
	ids = ['NC_045512.2', 'EPI_ISL_408669']
	parse_vapid_results(out_dir, ids)
	for i in ids:
		assert os.path.exists(f"{out_dir}/{i}/{i}.tsv")
		shutil.rmtree(f"{out_dir}/{i}/")