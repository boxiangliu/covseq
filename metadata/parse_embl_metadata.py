from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import click


def get_reference(record):
	for reference in record.annotations["references"]:
		if reference.title == ";":
			break
	return reference


def parse_assembly_data(assembly_data):
	container = {}
	for line in assembly_data:
		if line.startswith("Assembly Method"):
			container["Assembly_Method"] = line.split("::")[1].strip()
		elif line.startswith("Sequencing Technology"):
			container["Sequencing_Technology"] = line.split("::")[1].strip()
		elif line.startswith("Coverage"):
			container["Coverage"] = line.split("::")[1].strip()
	return container


def get_assembly_data_field(assembly_data, field):
	value = assembly_data[field] if field in assembly_data else ""
	return value


def get_qualifier_field(qualifiers, field):
	value = qualifiers[field][0] if field in qualifiers else ""
	return value


def parse_embl(embl):
	container = defaultdict(list)
	for record in embl:
		print(record.id)

		#####################
		# Submission fields #
		#####################
		reference = get_reference(record)
		assert reference.title == ";", \
			"Reference should be ';'!"
		
		container["Authors"].append(reference.authors)
		
		journal = reference.journal
		stuff, submitting_lab = journal.split(")", 1)
		container["Submitting_Lab"].append(submitting_lab.replace("to the INSDC.", "").strip())
		
		submission_date = stuff.split("(")[1].strip()
		container["Submission_Date"].append(submission_date)
		
		#####################
		# Sequencing fields #
		#####################
		if "comment" in record.annotations:
			assembly_data = record.annotations["comment"].split("\n")
			parsed_assembly_data = parse_assembly_data(assembly_data)
			for key in  ["Assembly_Method", "Sequencing_Technology", "Coverage"]:
				value = get_assembly_data_field(parsed_assembly_data, key)
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

	embl_df = pd.DataFrame(container)
	return embl_df


@click.command()
@click.option("--embl_fn", "-g", type=str, help="Input genbank file.")
@click.option("--out_fn", "-o", type=str, help="Output TSV file.")
def main(embl_fn, out_fn):
	print("#################")
	print("# EMBL metadata #")
	print("#################")
	print(f"EMBL: {embl_fn}")
	print(f"Output: {out_fn}")
	out_dir = os.path.dirname(out_fn)
	os.makedirs(out_dir, exist_ok=True)

	embl = SeqIO.parse(embl_fn, "embl")
	embl_df = parse_embl(embl)
	embl_df["Data_Source"] = "EMBL"
	embl_df.to_csv(out_fn, index=False, sep="\t")

if __name__ == "__main__":
	main()