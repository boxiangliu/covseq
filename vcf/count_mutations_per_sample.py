import glob
import sys
import os
import click
import pandas as pd
sys.path.append(".")
from utils import VCF, DefaultOrderedDict

def get_num_variants(vcf_dir):
	container = DefaultOrderedDict(list)
	for vcf_fn in glob.glob(f"{vcf_dir}/*.vcf.gz"):
		print(vcf_fn)
		ID = os.path.basename(vcf_fn).replace(".vcf.gz","")
		vcf = VCF(vcf_fn)
		num_variants = vcf.rowdata.shape[0]
		container["ID"].append(ID)
		container["num_variants"].append(num_variants)
	df = pd.DataFrame(container)
	return df

@click.command()
@click.option("--vcf_dir")
@click.option("--out_dir")
def main(vcf_dir, out_dir):
	os.makedirs(out_dir, exist_ok=True)
	num_variants = get_num_variants(vcf_dir)
	num_variants.to_csv(f"{out_dir}/mutations_per_sample.tsv", index=False, sep="\t")


if __name__ == "__main__":
	main()