import pandas as pd
import glob
import os
import click
from collections import defaultdict


def read_sample_info(in_dir):
	samples = []
	for fn in glob.glob(f"{in_dir}/*.info"):
		print(fn)
		sample = {}
		with open(fn) as f:
			for line in f:
				line = line.replace("\n", "")
				if "\t" not in line:
					continue
				split_line = line.split("\t")
				field = split_line[0].replace(":","")
				sample[field] = split_line[1]
		samples.append(sample)
	return samples


def samples2columns(samples):
	column_names = ["Accession ID",
	"Virus name",
	"Type",
	"Passage details/history",
	"Collection date",
	"Location",
	"Host",
	"Additional location information",
	"Gender",
	"Patient age",
	"Patient status",
	"Specimen source",
	"Additional host information",
	"Outbreak",
	"Last vaccinated",
	"Treatment",
	"Sequencing technology",
	"Assembly method",
	"Coverage",
	"Originating lab",
	"Sample ID given by the sample provider",
	"Submitting lab",
	"Sample ID given by the submitting laboratory",
	"Authors",
	"Submitter",
	"Submission Date"]

	columns = defaultdict(list)

	for s in samples:
		for c in column_names:
			if c in s:
				columns[c].append(s[c])
			else:
				columns[c].append("")

	return columns


@click.command()
@click.option("--metadata_dir", "-m", type=str, help="Directory where GISAID metadata is located.")
@click.option("--out_fn", "-o", type=str, help="Output file name.")
def main(metadata_dir, out_fn):
	# in_dir = "../data/gisaid/metadata/"
	# out_dir = "../processed_data/phenotype/"

	samples = read_sample_info(metadata_dir)
	columns = samples2columns(samples)
	df = pd.DataFrame(columns)
	df.to_csv(out_fn, sep="\t", index=False)

if __name__ == "__main__":
	main()