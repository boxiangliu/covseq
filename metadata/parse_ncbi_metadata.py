from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import click


def get_assembly_data_field(assembly_data, field):
	value = assembly_data[field] if field in assembly_data else ""
	return value


def get_qualifier_field(qualifiers, field):
	value = qualifiers[field][0] if field in qualifiers else ""
	return value


def parse_genbank(gb):
	container = defaultdict(list)
	for record in gb:
		print(record.id)

		#####################
		# Submission fields #
		#####################
		reference = record.annotations["references"][-1]
		assert reference.title == "Direct Submission", \
			"Reference should be direct submission!"
		
		container["Authors"].append(reference.authors)
		
		journal = reference.journal
		stuff, submitting_lab = journal.split(")", 1)
		container["Submitting_Lab"].append(submitting_lab.strip())
		
		submission_date = stuff.split("(")[1].strip()
		container["Submission_Date"].append(submission_date)
		
		#####################
		# Sequencing fields #
		#####################
		if "structured_comment" in record.annotations:
			assembly_data = record.annotations["structured_comment"]["Assembly-Data"]

			for field, key in [("Assembly Method", "Assembly_Method"), \
								("Sequencing Technology", "Sequencing_Technology"), \
								("Coverage", "Coverage")]:
				value = get_assembly_data_field(assembly_data, field)
				container[key].append(value)
		else:
			for field in ["Assembly_Method", "Sequencing_Technology", "Coverage"]:
				container[field].append("")

		################
		# Virus fields #
		################
		container["Accession_ID"].append(record.annotations["accessions"][0])
		container["Type"].append(record.annotations["organism"])
		container["Virus"].append(record.description)

		#################
		# Sample fields #
		#################
		source = record.features[0]
		assert source.type == "source", "Source type must be source!"
		qualifiers = source.qualifiers

		for field, key in [("isolation_source", "Specimen_Source"), \
						("host", "Host"), ("country", "Location"), \
						("collection_date", "Collection_Date")]:
			value = get_qualifier_field(qualifiers, field)
			container[key].append(value)

	genbank_df = pd.DataFrame(container)
	return genbank_df


def parse_csv(in_fn):
	df = pd.read_csv(in_fn)
	rename = {"Accession": "Accession_ID",
	"Release_Date": "Submission_Date",
	"Geo_Location": "Location", 
	"Isolation_Source": "Specimen_Source"}
	df.rename(columns=rename, inplace=True)
	df["Data_Source"] = "NCBI"
	df["Submission_Date"] = [x.split("T")[0] for x in df["Submission_Date"].tolist()]
	return df


@click.command()
@click.option("--gb_fn", "-g", type=str, help="Input genbank file.")
@click.option("--csv_fn", "-c", type=str, help="Input csv file.")
@click.option("--out_fn", "-o", type=str, help="Output TSV file.")
def main(gb_fn, csv_fn, out_fn):
	print("#################")
	print("# NCBI metadata #")
	print("#################")
	print(f"Genbank: {gb_fn}")
	print(f"CVS: {csv_fn}")
	print(f"Output: {out_fn}")
	out_dir = os.path.dirname(out_fn)
	os.makedirs(out_dir, exist_ok=True)

	if gb_fn:
		gb = SeqIO.parse(gb_fn, "genbank")
		genbank_df = parse_genbank(gb)
		genbank_df["Data_Source"] = "NCBI"
		genbank_df.to_csv(out_fn, index=False, sep="\t")
	elif csv_fn:
		csv_df = parse_csv(csv_fn)
		csv_df.to_csv(out_fn, index=False, sep="\t")


if __name__ == "__main__":
	main()