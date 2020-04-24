import pandas as pd
import sys
import click
import os
sys.path.append(".")
from utils import VCF, DefaultOrderedDict

def count_mutation(vcf):
	container = DefaultOrderedDict(list)
	for index, row in vcf.rowdata.iterrows():
		pos = row["POS"]
		ID = row["ID"]
		split_ID = ID.split(";")
		variant_count = len(split_ID)

		container["pos"].append(pos)
		container["variant_count"].append(variant_count)

	df = pd.DataFrame(container)
	return df

@click.command()
@click.option("--vcf_fn", help="VCF file.")
@click.option("--out_dir", help="Output directory.")
def main(vcf_fn, out_dir):
	os.makedirs(out_dir, exist_ok=True)
	vcf = VCF(vcf_fn)
	mutation_count = count_mutation(vcf)
	mutation_count.to_csv(f"{out_dir}/variant_count.tsv", index=False, sep="\t")

if __name__ == "__main__":
	main()