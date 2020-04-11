import os
import pytest
from fasta2vcf import *
from plot_mut_per_sample import *

def make_example_fastas():
	with open("ref.fasta", "w") as f:
		f.write(">ref\n")
		f.write("ATCGATCG\n")
		f.write("TCGAC\n")

	with open("qry.fasta", "w") as f:
		f.write(">qry\n")
		f.write("ATCAATCG\n")
		f.write("TCGAC\n")

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

def test_align():
	ref_fasta, qry_fasta = make_example_fastas()
	align(ref_fasta, qry_fasta, "test")

	assert os.path.exists("test.fasta")
	assert os.path.exists("test.ali")
	with open("test.ali") as f:
		line = f.readline()
		assert line.strip() == ">ref"
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


def test_align2vcf():
	ref = "ATCG"
	qry = "ATCG"
	rv, qv = align2vcf(ref, qry)
	assert rv == {}
	assert qv == {}

	ref = "AACG"
	qry = "ATCG"
	rv, qv = align2vcf(ref, qry)
	assert rv == {2:"A"}
	assert qv == {2:"T"}

	ref = "A-CG"
	qry = 'ATCG'
	rv, qv = align2vcf(ref, qry)
	assert rv == {1:"A-"}
	assert qv == {1:"AT"}

	ref = "ATCG"
	qry = "A-CG"
	rv, qv = align2vcf(ref, qry)
	assert rv == {1:"AT"}
	assert qv == {1:"A-"}


	ref = "ATTCG"
	qry = "A-CCG"
	rv, qv = align2vcf(ref, qry)
	assert rv == {1:"AT", 3:"T"}
	assert qv == {1:"A-", 3:"C"}


	ref = "ATT-G"
	qry = "A-CCG"
	rv, qv = align2vcf(ref, qry)
	assert rv == {1:"AT", 3:"T-"}
	assert qv == {1:"A-", 3:"CC"}


	ref = "ATT--"
	qry = "ATTCG"
	rv, qv = align2vcf(ref, qry)
	assert rv == {3:"T--"}
	assert qv == {3:"TCG"}


	ref = "--ATT"
	qry = "TCATT"
	rv, qv = align2vcf(ref, qry)
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


def test_get_mutation_count():
	count = get_mutation_count("../processed_data/fasta2vcf/EPI_ISL_402119.vcf")
	assert count == 1
	count = get_mutation_count("../processed_data/fasta2vcf/EPI_ISL_402127.vcf")
	assert count == 3
