import sys
sys.path.append(sys.path.abspath("."))
from utils import VCF, DefaultOrderedDict
from random import seed, randint
import numpy as np
from scipy.stats import pearsonr
import plotly.express as px
import os
import pickle
import pandas as pd

vcf_fn = "../data/aggregated/vcf/merged/filtered.vcf.gz"
out_dir = "../processed_data/mutation_distribution/mutation_correlation/"
os.makedirs(out_dir, exist_ok=True)

def calc_r(snp1, snp2):
	return pearsonr(snp1, snp2)[0]

def standardize(x):
	mu = np.mean(x)
	sd = np.std(x)
	return (x - mu)/sd

print("Getting VCF.")
if not os.path.exists(f"{out_dir}/vcf.pkl"):
	vcf = VCF(vcf_fn)
	with open(f"{out_dir}/vcf.pkl", "wb") as f:
		pickle.dump(vcf, f)
else:
	with open(f"{out_dir}/vcf.pkl", "rb") as f:
		vcf = pickle.load(f)


n_variant = vcf.rowdata.shape[0]
n_sample = vcf.data.shape[1]
vcf.data.index = vcf.rowdata["POS"]

pos = vcf.rowdata.loc[:,"POS"]
genotype = vcf.data.to_numpy()
allele_freq = vcf.data.apply(np.mean, axis=1)

print("Getting R2.")
genotype_std = np.apply_along_axis(standardize, 1, genotype)
r = genotype_std @ genotype_std.T / n_sample
r2 = r**2


if not os.path.exists(f"{out_dir}/r2_long.pkl"):
	container = DefaultOrderedDict(list)
	for i in range(n_variant):
		print(f"First SNP: {i}")
		pos_i = pos.iloc[i]
		for j in range(i+1, n_variant):
			print(f"Second SNP: {j}")
			pos_j = int(pos.iloc[j])
			dist_ij = pos_j - pos_i
			container["pos_1"].append(pos_i)
			container["pos_2"].append(pos_j)
			container["dist"].append(dist_ij)
			container["r2"].append(r2[i,j])

	r2_long = pd.DataFrame(container)

	with open(f"{out_dir}/r2_long.pkl", "wb") as f
		pickle.dump(r2_long, f)
else:
	r2_long = pickle.load(f"{out_dir}/r2_long.pkl")

threshold = 0.001
r2_long["r2"].describe(percentiles=np.arange(0.1,1,0.01))
r2_filt = r2_long[r2_long["r2"] > threshold]

fig_scat = px.scatter(r2_filt, x="dist", y="r2")
fig_scat.show()
fig_hist = px.histogram(r2_filt, x="r2")
fig_hist.show()
r2_long[(r2_long["pos_1"] == 8782) & (r2_long["pos_2"] == 28144)]

# Goal: Look at variant pairs with very high correlation (r2 = 1)
print("Investigate variant pairs with r2=1.")
r2_eq_1 = r2_long[r2_long["r2"] == 1]

r2_eq_1
#           pos_1  pos_2   dist   r2
# 13062253   4706  14318   9612  1.0
# 23349644  10533  29460  18927  1.0
# 37230154  23045  27786   4741  1.0
# 37463825  23485  24087    602  1.0
# 37464437  23485  26005   2520  1.0
# 37588226  23665  24769   1104  1.0
# 37862636  24087  26005   1918  1.0

# Goal: Look at a few example variants.
allele_freq.loc[23665]
# 2.8544515171409814e-05
allele_freq.loc[24769]
# 2.8544515171409814e-05
allele_freq.describe()
# count    8908.000000
# mean        0.000530
# std         0.012836
# min         0.000029
# 25%         0.000029
# 50%         0.000057
# 75%         0.000114
# max         0.675620
# dtype: float64
# Conclusion: The allele frequency of these 
# variants are very low. In fact, these variants 
# have the lowest allele frequency.

# Question: Are the variants singletons?
1/n_sample 
# 2.8544515171409814e-05
# Conclusion: Yes! These variants are singletons!

# Questions: Which strain carry these variants? 
EPI_ISL_420357 = vcf.data.loc[:,vcf.data.loc[4706] == 1].iloc[:,0]
EPI_ISL_420357.loc[EPI_ISL_420357==1]
# POS
# 4706     1
# 14318    1
# 23403    1
# 24389    1
allele_freq.loc[23403]
# 0.6756201295920988
allele_freq.loc[24389]
# 0.0006279793337710159
# Name: EPI_ISL_420357, dtype: int64
# Conclusion: Variants 4706 and 14318 are singletons 
# that only exists in EPI_ISL_420357
# The other two variants 23403 and 24389 have higher
# allele frequencies. They are probably ancestral 
# alleles passed onto EPI_ISL_420357.

# Goal: Remove singletons.
def remove_singletons(vcf):
	n_alt = vcf.data.apply(sum, axis=1)
	vcf.data = vcf.data.

# seed(42)
# n = 1000
# snp1 = np.array([randint(0,1) for x in range(n)])
# snp2 = np.array([randint(0,1) for x in range(n)])
# snp1_std = standardize(snp1)
# snp2_std = standardize(snp2)
# r = snp1_std @ snp2_std / n






