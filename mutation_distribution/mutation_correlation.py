import sys
sys.path.append(".")
from utils import VCF
from random import seed, randint
import numpy as np
from collections import OrderedDict
from scipy.stats import pearsonr
import plotly.express as px

vcf_fn = "../data/aggregated/vcf/merged/filtered.vcf.gz"
vcf = VCF(vcf_fn)

n_variant = vcf.rowdata.shape[0]
n_sample = vcf.data.shape[1]

pos = vcf.rowdata.loc[:,"POS"]
genotype = vcf.data.to_numpy()

# container = OrderedDict()
# for i in range(nrow):
# 	print(f"First SNP: {i}")
# 	pos_i = pos.iloc[i]
# 	snp_i = vcf.data.iloc[i,:]
# 	for j in range(i, nrow):
# 		print(f"Second SNP: {j}")
# 		pos_j = int(pos.iloc[j])
# 		snp_j = vcf.data.iloc[j,:]
# 		container[(pos_i, pos_j)] = calc_r(snp_i, snp_j)


def calc_r(snp1, snp2):
	return pearsonr(snp1, snp2)[1]


def standardize(x):
	mu = np.mean(x)
	sd = np.std(x)
	return (x - mu)/sd


genotype_std = np.apply_along_axis(standardize, 1, genotype)
r = genotype_std @ genotype_std.T / n_sample
fig = px.imshow(r)
fig.show()


# seed(42)
# n = 1000
# snp1 = np.array([randint(0,1) for x in range(n)])
# snp2 = np.array([randint(0,1) for x in range(n)])
# snp1_std = standardize(snp1)
# snp2_std = standardize(snp2)
# r = snp1_std @ snp2_std / n






