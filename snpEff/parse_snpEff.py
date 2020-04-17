import sys
sys.path.append(".")
from utils import VCF

in_fn = "../data/aggregated/vcf/merged/annotated.vcf.gz"
vcf = VCF(in_fn)