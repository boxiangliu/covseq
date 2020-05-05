import os
import click
import pandas as pd
import sys
sys.path.append(".")
from utils import VCF, DefaultOrderedDict


HEADERS = ["ALLELE", "EFFECT", "IMPACT", \
	"GENE", "GENEID", "FEATURE", "FEATUREID", \
	"BIOTYPE", "RANK", "HGVS_C", "HGVS_P", \
	"CDNA_POS", "CDS_POS", \
	"AA_POS", "DISTANCE", "ERRORS"]

def parse_snpEff(vcf_fn, ORF1a=None):
	"""
	vcf_fn: VCF file
	ORF1a: a list of ORF1a start and end, i.e. [266, 13483]. Providing this parameter will 
		force the GENE column to append ",ORF1a" if the variant falls within ORF1a. 
	"""
	vcf = VCF(vcf_fn)
	container = DefaultOrderedDict(list)
	container["POS"] = vcf.rowdata["POS"].tolist()
	container["REF"] = vcf.rowdata["REF"].tolist()
	container["ALT"] = vcf.rowdata["ALT"].tolist()

	ANN_0 = [x.split(",")[0].replace("ANN=","") for x in vcf.rowdata["INFO"]]
	for ann in ANN_0:
		split_ann = ann.split("|")
		assert len(split_ann) == 16
		for i, field in enumerate(split_ann):
			container[HEADERS[i]].append(field)
	assert container["ALT"] == container["ALLELE"]
	del container["ALLELE"]
	
	if ORF1a:
		ORF1a_start = ORF1a[0]
		ORF1a_end = ORF1a[1]
		for i, pos in enumerate(container["POS"]):
			if ORF1a_start <= pos and pos <= ORF1a_end:
				container["GENE"][i] += ",ORF1a"

	parsed = pd.DataFrame(container)
	return parsed


@click.command()
@click.option("--vcf_fn", "-v", type=str, help="VCF file input.")
@click.option("--out_fn", "-o", type=str, help="Output file.")
def main(vcf_fn, out_fn):
	snpEff = parse_snpEff(vcf_fn)
	out_dir = os.path.dirname(out_fn)
	os.makedirs(out_dir, exist_ok=True)
	snpEff.to_csv(out_fn, index=False, sep="\t")

if __name__ == "__main__":
	main()