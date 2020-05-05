import os
import pandas as pd
import click

in_fn = "../data/cngb/metadata/CNGBdb_VirusDIP_excel20200411_all(24)_57b4ce53c4d6c49c978596677a112211.csv"
out_fn = "../data/aggregated/metadata/individual/cngb.tsv"

@click.command()
@click.option("--in_fn", "-i", type=str, help="Input CNGB file.")
@click.option("--out_fn", "-o", type=str, help="Output CNGB file.")
def main(in_fn, out_fn):
	print("#################")
	print("# CNGB metadata #")
	print("#################")
	print(f"CNGB: {in_fn}")
	print(f"Output: {out_fn}")
	out_dir = os.path.dirname(out_fn)
	os.makedirs(out_dir, exist_ok=True)

	df = pd.read_csv(in_fn)
	old2new = {"Sequence ID": "Accession_ID",
		"Organism": "Type",
		"Sample collection date": "Collection_Date",
		"Released date": "Submission_Date", 
		"Data source platform": "Data_Source",
		"Submitter organization": "Submitting_Lab", 
		"Sequencing technology/Platform": "Sequencing_Technology",
		"Assembly method": "Assembly_Method", 
		"Virus name": "Virus"}

	df.rename(columns=old2new, inplace=True)
	df.drop(["Tax ID","Length","Originating Lab", "Literature", "Files", "Browse"], axis=1, inplace=True)
	df.to_csv(out_fn, sep="\t", index=False)

if __name__ == "__main__":
	main()
