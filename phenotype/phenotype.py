import pandas as pd
import glob
import os
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



def main():
	in_dir = "../data/gisaid/"
	out_dir = "../processed_data/phenotype/"
	os.makedirs(out_dir, exist_ok=True)

	samples = read_sample_info(in_dir)
	columns = samples2columns(samples)
	df = pd.DataFrame(columns)
	df.to_csv(f"{out_dir}/phenotype.tsv", sep="\t", index=False)

if __name__ == "__main__":
	main()