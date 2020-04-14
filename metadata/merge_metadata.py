import pandas as pd
import click

SOURCES = ["gisaid","ncbi","embl","cngb"]
COLUMNS = ["Accession_ID",
	"Virus",
	"Type",
	"Data_Source",
	"Collection_Date",
	"Location",
	"Host",
	"Patient_Gender",
	"Patient_Age",
	"Patient_Status",
	"Specimen_Source",
	"Sequencing_Technology",
	"Assembly_Method",
	"Coverage",
	"Submitting_Lab",
	"Authors",
	"Submitter",
	"Submission_Date"]

@click.command()
@click.option("--in_dir", "-i", type=str, help="Input directory.")
@click.option("--out_fn", "-o", type=str, help="Output file.")
def main(in_dir, out_fn):
	container = []
	for s in SOURCES:
		fn = f"{in_dir}/{s}.tsv"
		df = pd.read_table(fn)
		for c in COLUMNS:
			if c not in df.columns:
				df[c] = ""
		df = df[COLUMNS]
		container.append(df)

	concat = pd.concat(container)
	concat.to_csv(out_fn, index=False, sep="\t")

if __name__ == "__main__":
	main()