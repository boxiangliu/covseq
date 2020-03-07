"""
Annotate virus sequence using VAPID
"""
import os
from Bio import SeqIO
import pytest
from datetime import date
import subprocess
from parse_vapid import Annotation

def get_fasta_ids(fasta_fn):
	records = SeqIO.parse(fasta_fn, format="fasta")
	ids = []
	for record in records:
		ids.append(record.id)
	return ids


def create_metadata(ids, out_dir):
	today = str(date.today())
	region = "Unknown"
	coverage = str(0)

	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	with open(f"{out_dir}/metadata.csv", "w") as f:
		f.write("strain,collection-date,country,coverage\n")

		for id_ in ids:
			out = ",".join([id_, today, region, coverage])
			f.write(out + "\n")


def mafft_is_installed():
	try:
		subprocess.run(["mafft", "--version"]) 
		return True
	except FileNotFoundError:
		return False


def tbl2asn_is_installed():
	try:
		subprocess.run(["tbl2asn", "--version"])
		return True
	except:
		return False


def call_vapid(fasta_fn, metadata_fn, out_dir):
	if not mafft_is_installed():
		raise Exception("mafft is not installed. Visit https://mafft.cbrc.jp/alignment/software/")
	
	if not tbl2asn_is_installed():
		raise Exception("tbl2asn is not installed. Visit https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/")

	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	ids = get_fasta_ids(fasta_fn)

	cmd = f"python3 src/vapid_minimal/vapid3.py \
	--db src/vapid_minimal/all_virus.fasta \
	--r NC_045512.2 \
	--metadata_loc {metadata_fn} \
	{fasta_fn} \
	src/vapid_minimal/template.sbt"

	cmd = cmd.split()
	output = subprocess.run(cmd, capture_output=True)

	cwd = os.getcwd()
	for i in ids:
		os.rename(f"{cwd}/{i}", f"{cwd}/{out_dir}/{i}")

	return output


def parse_vapid_results(out_dir, ids):
	for i in ids:
		prefix = f"{out_dir}/{i}/{i}"
		anno = Annotation(prefix)
		anno.anno_df.to_csv(prefix + ".tsv", sep="\t", index=False)