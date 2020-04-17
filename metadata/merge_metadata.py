import pandas as pd
import click
import dateutil.parser
import re
import sys
sys.path.append(".")
from utils import VCF

SOURCES = ["gisaid_acknowledgement","ncbi","embl","cngb"]
COLUMNS = ["Accession_ID",
	"Virus",
	"Data_Source",
	"Collection_Date",
	"Location",
	"Submitting_Lab",
	"Authors"]


def format_dates(dates):
	formatted = []
	for d in dates:
		try:
			parsed_date = dateutil.parser.parse(d)
			split_date = re.split(" |-", d)
			if len(split_date) == 3:
				formatted.append(parsed_date.strftime("%Y-%m-%d"))
			elif len(split_date) == 2:
				formatted.append(parsed_date.strftime("%Y-%m"))
			elif len(split_date) == 1:
				formatted.append(parsed_date.strftime("%Y"))
			else:
				raise Exception("Please check the date.")
		except TypeError:
			formatted.append("")
	return formatted


def parse_location(locations):
	countries = []
	regions = []
	for l in locations:
		try: 
			l = l.replace("Europe /", "", 1).\
			replace("Europe/", "", 1).\
			replace("Asia /", "", 1).\
			replace("Asia/", "", 1).\
			replace("North America /", "", 1).\
			replace("North America/", "", 1).\
			replace("Oceania /", "", 1).\
			replace("Oceania/", "", 1).\
			replace("Africa /", "", 1).\
			replace("Africa/", "", 1).\
			replace("South America /", "", 1).\
			replace("South America/", "", 1).\
			replace("Central America /", "", 1).\
			replace("Central America/", "", 1)
			split_location = re.split(":|/", l)
			country = split_location[0].strip()
			if len(split_location) == 0:
				country = ""
				region = ""
			elif len(split_location) == 1:
				country = split_location[0].strip()
				region = ""
			else:
				country = split_location[0].strip()
				region = split_location[1].strip()

			if country == "Hong Kong":
				country = "China"
				region = "Hong Kong"

			if country == "England":
				country = "United Kingdom"
				region = "England"

		except AttributeError:
			country = ""
			region = ""

		countries.append(country)
		regions.append(region)
	return countries, regions


@click.command()
@click.option("--in_dir", "-i", type=str, help="Input directory.")
@click.option("--out_prefix", "-o", type=str, help="Output file prefix.")
@click.option("--vcf_fn", "-v", type=str, help="A VCF file. If set, the program will output an additional file with only entries that appear in the VCF.")
def main(in_dir, out_prefix, vcf_fn):
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

	# Format dates
	dates = concat["Collection_Date"].tolist()
	formatted = format_dates(dates)
	concat["Collection_Date"] = formatted

	# Format locations.
	location = concat["Location"].tolist()
	country, region = parse_location(location)
	concat["Country"] = country
	concat["Region"] = region
	concat.to_csv(f"{out_prefix}.tsv", index=False, sep="\t")

	if vcf_fn:
		vcf = VCF(vcf_fn)
		concat_in_vcf = concat[concat["Accession_ID"].isin(vcf.data.columns)]
		concat_in_vcf = concat_in_vcf.loc[~concat_in_vcf["Accession_ID"].duplicated(),:]
		concat_in_vcf.to_csv(f"{out_prefix}_in_vcf.tsv", index=False, sep="\t")

if __name__ == "__main__":
	main()


