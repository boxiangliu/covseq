import pytest
import shutil
import os
from annotate import *

def test_get_fasta_ids():
	ids = get_fasta_ids("example.fasta")
	assert ids == ['NC_045512.2', 'EPI_ISL_408669']

def test_create_metadata():
	metadata_fn = create_metadata(["A","B"], out_dir="results/example/")
	assert os.path.normpath(metadata_fn) == "results/example/metadata.csv"
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
	parse_vapid(out_dir, ids)
	for i in ids:
		assert os.path.exists(f"{out_dir}/{i}/{i}.tsv")


def test_clean():
	in_dir = "results/example/vapid/"
	out_dir = "results/example/annotation/"
	clean(in_dir, out_dir)
	assert not os.path.exists(in_dir)
	assert os.path.exists(f"{out_dir}/NC_045512.2/NC_045512.2.tsv") 
	assert os.path.exists(f"{out_dir}/EPI_ISL_408669/EPI_ISL_408669.tsv") 
	shutil.rmtree(out_dir)


def test_annotate():
	fasta_fn = "example.fasta"
	out_dir = "results/example/"
	ids = ['NC_045512.2', 'EPI_ISL_408669']
	annotate(fasta_fn, out_dir)
	for i in ids:
		assert os.path.exists(f"{out_dir}/annotation/{i}/{i}.tsv")
	shutil.rmtree(out_dir)

