import os
from collections import OrderedDict
import pytest
import subprocess
from fasta2vcf import *
from filter_samples import *

def make_example_fastas():
	with open("ref.fasta", "w") as f:
		f.write(">1\n")
		f.write("ATCGATCG\n")
		f.write("TCGAC\n")

	with open("qry.fasta", "w") as f:
		f.write(">qry\n")
		f.write("ATCAATCG\n")
		f.write("TCGAC\n")

	cmd = "samtools faidx ref.fasta"
	subprocess.run(cmd, shell=True)

	return "ref.fasta", "qry.fasta"


def remove_example_fastas():
	os.remove("qry.fasta")
	os.remove("ref.fasta")


def make_example_align():
	with open("test.ali", "w") as f:
		f.write(">ref\n")
		f.write("ATCGATCG\n")
		f.write("TCGAC\n")
		f.write(">qry\n")
		f.write("ATCAATCG\n")
		f.write("TCGAC\n")

	return "test.ali"

def test_get_IDs_from_align_file():
	make_example_align()
	ref_id, qry_id = get_IDs_from_align_file("test.ali")
	assert ref_id == "ref"
	assert qry_id == "qry"


def test_align():
	ref_fasta, qry_fasta = make_example_fastas()
	align(ref_fasta, qry_fasta, "test")

	assert os.path.exists("test.fasta")
	assert os.path.exists("test.ali")
	with open("test.ali") as f:
		line = f.readline()
		assert line.strip() == ">1"
		line = f.readline()
		assert line.strip() == "ATCGATCGTCGAC".lower()
		line = f.readline()
		assert line.strip() == ">qry"
		line = f.readline()
		assert line.strip() == "ATCAATCGTCGAC".lower()

	remove_example_fastas()
	os.remove("test.ali")
	os.remove("test.fasta")


def test_parse_align():
	align_fn = make_example_align()

	ref, qry = parse_align(align_fn)
	assert ref == "ATCGATCGTCGAC"
	assert qry == "ATCAATCGTCGAC"

	os.remove("test.ali")


def test_align2variant():
	ref = "ATCG"
	qry = "ATCG"
	rv, qv = align2variant(ref, qry)
	assert rv == {}
	assert qv == {}

	ref = "AACG"
	qry = "ATCG"
	rv, qv = align2variant(ref, qry)
	assert rv == {2:"A"}
	assert qv == {2:"T"}

	ref = "A-CG"
	qry = 'ATCG'
	rv, qv = align2variant(ref, qry)
	assert rv == {1:"A-"}
	assert qv == {1:"AT"}

	ref = "ATCG"
	qry = "A-CG"
	rv, qv = align2variant(ref, qry)
	assert rv == {1:"AT"}
	assert qv == {1:"A-"}


	ref = "ATTCG"
	qry = "A-CCG"
	rv, qv = align2variant(ref, qry)
	assert rv == {1:"AT", 3:"T"}
	assert qv == {1:"A-", 3:"C"}


	ref = "ATT-G"
	qry = "A-CCG"
	rv, qv = align2variant(ref, qry)
	assert rv == {1:"AT", 3:"T-"}
	assert qv == {1:"A-", 3:"CC"}


	ref = "ATT--"
	qry = "ATTCG"
	rv, qv = align2variant(ref, qry)
	assert rv == {3:"T--"}
	assert qv == {3:"TCG"}


	ref = "--ATT"
	qry = "TCATT"
	rv, qv = align2variant(ref, qry)
	assert rv == {0:"--"}
	assert qv == {0:"TC"}

def test_save_vcf():
	ref_variant = OrderedDict({1:"A", 3:"T"})
	qry_variant = OrderedDict({1:"C", 3:"TG"})
	qry_name = "test"
	out_fn = "test.vcf"
	save_vcf(ref_variant, qry_variant, qry_name, out_fn)
	assert os.path.exists("test.vcf")
	os.remove("test.vcf")


def test_align2vcf_and_get_mutation_count():
	ref_fasta, qry_fasta = make_example_fastas()
	align(ref_fasta, qry_fasta, "test")
	align2vcf("test.ali", ref_fasta, qry_id = "qry", out_prefix="test", compress_vcf=True, clean_up=True, verbose=True)
	assert os.path.exists("test.vcf.gz")

	count = get_mutation_count("test.vcf.gz")
	assert count == 1

	align(ref_fasta, qry_fasta, "test")
	align2vcf("test.ali", ref_fasta, qry_id = "qry", out_prefix="test", compress_vcf=False, clean_up=True)
	assert os.path.exists("test.vcf")

	count = get_mutation_count("test.vcf")
	assert count == 1

	os.remove("qry.fasta")
	os.remove("ref.fasta")
	os.remove("ref.fasta.fai")
	os.remove("test.vcf.gz")
	os.remove("test.vcf.gz.tbi")
	os.remove("test.vcf")