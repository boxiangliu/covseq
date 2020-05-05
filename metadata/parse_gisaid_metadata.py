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
	column_names = {"Accession ID": "Accession_ID",
	"Virus name": "Virus",
	"Type": "Type",
	"Passage details/history": "Passage_Details/History",
	"Collection date": "Collection_Date",
	"Location": "Location",
	"Host": "Host",
	"Additional location information": "Additional_Location_Information",
	"Gender": "Patient_Gender",
	"Patient age": "Patient_Age",
	"Patient status": "Patient_Status",
	"Specimen source": "Specimen_Source",
	"Additional host information": "Additional_Host_Information",
	"Outbreak": "Outbreak",
	"Last vaccinated": "Last_Vaccinated",
	"Treatment": "Treatment",
	"Sequencing technology": "Sequencing_Technology",
	"Assembly method": "Assembly_Method",
	"Coverage": "Coverage",
	"Originating lab": "Originating_Lab",
	"Sample ID given by the sample provider": "Sample_ID_Given_by_the_Sample_Provider",
	"Submitting lab": "Submitting_Lab",
	"Sample ID given by the submitting laboratory": "Sample_ID_Given_by_the_Submitting_Laboratory",
	"Authors": "Authors",
	"Submitter": "Submitter",
	"Submission Date": "Submission_Date"}

	columns = defaultdict(list)

	for s in samples:
		for k, v in column_names.items():
			if k in s:
				columns[v].append(s[k])
			else:
				columns[v].append("")

	return columns


@click.command()
@click.option("--metadata_dir", "-m", type=str, help="Directory where GISAID metadata is located.")
@click.option("--acknowledgement_fn", "-m", type=str, help="GISAID acknowledgement file.")
@click.option("--out_fn", "-o", type=str, help="Output file name.")
@click.option("--type", "-t", type=str, help="Type of metadata file (detail, acknowledgement).")
def main(metadata_dir, acknowledgement_fn, out_fn, type):
	print("###################")
	print("# GISAID metadata #")
	print("###################")
	print(f"GISAID: {metadata_dir}")
	print(f"Output: {out_fn}")

	if type == "detail":
		samples = read_sample_info(metadata_dir)
		columns = samples2columns(samples)
		df = pd.DataFrame(columns)
		df["Data_Source"] = "GISAID"
		df.to_csv(out_fn, sep="\t", index=False)
	elif type == "acknowledgement":
		# metadata_dir = "../data/gisaid/metadata/acknowledgement/"
		fn = acknowledgement_fn
		df = pd.read_excel(fn, skiprows=[0,1,3])
		df["Data_Source"] = "GISAID"
		df.rename(columns={"Accession ID": "Accession_ID", \
			"Virus name": "Virus", "Collection date": "Collection_Date", \
			"Originating lab": "Originating_Lab", \
			"Submitting lab": "Submitting_Lab"}, inplace=True)
		df.to_csv(out_fn, sep="\t", index=False)
	else:
		raise Exception(f"{type} is not recognized.")

if __name__ == "__main__":
	main()