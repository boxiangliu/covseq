import datetime
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
sys.path.append(".")
from utils import VCF

CDS = {"5'UTR": [1,265],
"ORF1ab": [266, 21555],
"S": [21563, 25384], 
"ORF3a": [25393, 26220], 
"E": [26245, 26472], 
"M": [26523, 27191], 
"ORF6": [27202, 27387], 
"ORF7a": [27394, 27759],
"ORF7b": [27756, 27887], 
"ORF8": [27894, 28259],
"N": [28274, 29533],
"ORF10": [29558, 29674],
"3'UTR": [29675, 29903]}


def get_CDS_length(CDS):
	CDS_length = {}
	total_length = 29903
	for k, (start, end) in CDS.items():
		CDS_length[k] = end - start
	CDS_length["intergenic"] = total_length - sum(CDS_length.values())
	return CDS_length


def validate(date):
	try:
		datetime.datetime.strptime(date,"%Y-%m-%d")
		return True
	except:
		return False


def get_sorted_date(vcf):
	# human_data = vcf.coldata[vcf.coldata["Host"] == "Human"]
	human_data = vcf.coldata
	valid_date = [validate(x) for x in human_data["Collection_Date"].tolist()]
	sample_date = human_data[valid_date][["Accession_ID", "Collection_Date"]]
	sample_date.sort_values("Collection_Date", inplace=True)
	return sample_date


def get_variant_CDS(vcf, CDS):
	variant_CDS = {}
	for _, row in vcf.rowdata.iterrows():
		ref = row["REF"]
		pos = row["POS"]
		alt = row["ALT"]
		key = f'{ref}{pos}{alt}'
		for c, (start, end) in CDS.items():
			if start <= pos and pos <= end:
				variant_CDS[key] = c
				break 
		if key not in variant_CDS:
			variant_CDS[key] = "intergenic"
	return variant_CDS


def get_mutation_first_apparence(vcf, sample_date, variant_CDS):
	mutation_date = {}
	for _, s_row in sample_date.iterrows():
		ID = s_row["Accession_ID"]
		date = s_row["Collection_Date"]
		genotypes = vcf.data.loc[:,ID]
		variants = vcf.rowdata[genotypes.isin([1])]
		for _, v_row in variants.iterrows():
			key = f'{v_row["REF"]}{v_row["POS"]}{v_row["ALT"]}'
			if key not in mutation_date:
				mutation_date[key] = (date, variant_CDS[key])

	mutations = []
	dates = []
	cds = []
	for k, (d, c) in mutation_date.items():
		mutations.append(k)
		dates.append(d)
		cds.append(c)
	
	mutation_date_cds = pd.DataFrame({"mutation": mutations, "date": dates, "cds": cds})
	return mutation_date_cds


def get_new_mutations_each_day(mutation_date_cds):
	mutation_per_day = mutation_date_cds.groupby("date").agg(count=pd.NamedAgg(column="mutation", aggfunc="count")).reset_index()
	mutation_per_day["cumsum"] = mutation_per_day["count"].cumsum()
	return mutation_per_day


def get_new_mutations_each_day_strat_by_CDS(mutation_date_cds, CDS_length):
	mutation_per_day_strat_by_CDS = mutation_date_cds.groupby(["date", "cds"]).agg(count=pd.NamedAgg(column="mutation", aggfunc="count"))
	mutation_per_day_strat_by_CDS["cumsum"] = mutation_per_day_strat_by_CDS.groupby("cds").cumsum()
	mutation_per_day_strat_by_CDS = mutation_per_day_strat_by_CDS.reset_index()
	mutation_per_day_strat_by_CDS["count_per_length"] = mutation_per_day_strat_by_CDS.apply(lambda x: x["count"]/CDS_length[x["cds"]], axis=1)
	mutation_per_day_strat_by_CDS["cumsum_per_length"] = mutation_per_day_strat_by_CDS.apply(lambda x: x["cumsum"]/CDS_length[x["cds"]], axis=1)
	return mutation_per_day_strat_by_CDS


def fill_in_missing_days(mutation_per_day_strat_by_CDS):
	mutation_per_day_strat_by_CDS_wide = mutation_per_day_strat_by_CDS.pivot(index="date", columns="cds", values="count_per_length").fillna(0)
	df1 = mutation_per_day_strat_by_CDS_wide.reset_index().melt(id_vars="date", var_name="cds", value_name="count_per_length")
	
	proto_df2 = {}
	for index, column in mutation_per_day_strat_by_CDS_wide.iteritems():
		proto_df2[index] = column.cumsum()

	df2 = pd.DataFrame(proto_df2).reset_index()
	df2 = df2.melt(id_vars="date", var_name="cds", value_name="cumsum_per_length")
	df3 = pd.merge(df1, df2, on=["date","cds"])

	return df3


def plot_new_mutations_over_time(mutation_per_day, out_fn):
	fig, ax = plt.subplots(figsize=(16,4))
	font = {'family': 'sans-serif', 'color':  'black', 'weight': 'bold', 'size': 12}
	ax = sns.scatterplot(data=mutation_per_day, x="date", y="count", ax=ax)
	ax.set_xticklabels(labels=mutation_per_day["date"],rotation=45, ha='right')
	ax.set_xlabel(None)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_ylabel("Number of Novel Mutations", fontdict=font)
	fig.tight_layout()
	fig.savefig(out_fn)


def plot_cumulative_mutations_over_time(mutation_per_day, out_fn):
	fig, ax = plt.subplots(figsize=(16,4))
	font = {'family': 'sans-serif', 'color':  'black', 'weight': 'bold', 'size': 12}
	ax = sns.scatterplot(data=mutation_per_day, x="date", y="cumsum", ax=ax)
	ax.set_xticklabels(labels=mutation_per_day["date"],rotation=45, ha='right')
	ax.set_xlabel(None)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_ylabel("Cum. Number of Novel Mutations", fontdict=font)
	fig.tight_layout()
	fig.savefig(out_fn)


def plot_new_mutations_per_CDS_over_time(mutation_per_day, out_fn):
	fig, ax = plt.subplots(figsize=(16,4))
	font = {'family': 'sans-serif', 'color':  'black', 'weight': 'bold', 'size': 12}
	ax = sns.lineplot(data=mutation_per_day, x="date", y="count_per_length", hue="cds")
	ax.set_xticklabels(labels=mutation_per_day["date"].unique(),rotation=45, ha='right')
	ax.set_xlabel(None)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_ylabel("Number of Novel Mutations\nper Nuceotide", fontdict=font)
	fig.tight_layout()
	fig.savefig(out_fn)


def plot_cumulative_mutations_per_CDS_over_time(mutation_per_day_strat_by_CDS, out_fn):
	fig, ax = plt.subplots(figsize=(16,4))
	font = {'family': 'sans-serif', 'color':  'black', 'weight': 'bold', 'size': 12}
	ax = sns.lineplot(data=mutation_per_day_strat_by_CDS, x="date", y="cumsum_per_length", hue="cds")
	ax.set_xticklabels(labels=mutation_per_day_strat_by_CDS["date"].unique(),rotation=45, ha='right')
	ax.set_xlabel(None)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_ylabel("Cum. Number of Novel \nMutations per Nuceotide", fontdict=font)
	fig.tight_layout()
	fig.savefig(out_fn)


def plot_cumulative_mutations_per_CDS(mutation_per_day_strat_by_CDS, out_fn):
	latest_date = mutation_per_day_strat_by_CDS["date"].max()
	to_plot = mutation_per_day_strat_by_CDS[mutation_per_day_strat_by_CDS["date"]==latest_date]
	to_plot.sort_values("cumsum_per_length", inplace=True)
	fig, ax = plt.subplots(figsize=(16,4))
	font = {'family': 'sans-serif', 'color':  'black', 'weight': 'bold', 'size': 12}
	ax = sns.barplot(data=to_plot, x="cds", y="cumsum_per_length")
	ax.set_xlabel(None)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_ylabel("Cum. Number of Novel \nMutations per Nuceotide", fontdict=font)
	fig.tight_layout()
	fig.savefig(out_fn)


def main():
	vcf_fn = "../data/aggregated/vcf/merged/filtered.vcf.gz"
	pheno_fn = "../data/aggregated/metadata/merged.tsv"
	# pheno_fn = "../processed_data/phenotype/phenotype.tsv"
	out_dir = "../processed_data/mutation_distribution/mutation_over_time/"
	os.makedirs(out_dir, exist_ok=True)

	print("Reading VCF.")
	vcf = VCF(vcf_fn, pheno_fn)

	print("Getting CDS.")
	variant_CDS = get_variant_CDS(vcf, CDS)
	CDS_length = get_CDS_length(CDS)

	print("Getting date.")
	sample_date = get_sorted_date(vcf)
	mutation_date_cds = get_mutation_first_apparence(vcf, sample_date, variant_CDS)
	mutation_per_day = get_new_mutations_each_day(mutation_date_cds)
	mutation_per_day_strat_by_CDS = get_new_mutations_each_day_strat_by_CDS(mutation_date_cds, CDS_length)
	mutation_per_day_strat_by_CDS = fill_in_missing_days(mutation_per_day_strat_by_CDS)

	print("Plotting.")
	plot_new_mutations_over_time(mutation_per_day, f"{out_dir}/mutation_count_per_day.pdf")
	plot_cumulative_mutations_over_time(mutation_per_day, f"{out_dir}/mutation_cumsum_per_day.pdf")

	plot_new_mutations_per_CDS_over_time(mutation_per_day_strat_by_CDS, f"{out_dir}/mutation_count_per_CDS_per_day.pdf")
	plot_cumulative_mutations_per_CDS_over_time(mutation_per_day_strat_by_CDS, f"{out_dir}/mutation_sumsum_per_CDS_per_day.pdf")
	plot_cumulative_mutations_per_CDS(mutation_per_day_strat_by_CDS,f"{out_dir}/mutation_sumsum_per_CDS.pdf")

# to_plot = mutation_per_day_strat_by_CDS[mutation_per_day_strat_by_CDS["cds"]=="ORF1ab"]
# plt.plot(to_plot["date"], to_plot["cumsum"])
# locs, labels = plt.xticks()
# plt.xticks(locs, labels, rotation=45, ha='right')

# to_plot = mutation_per_day_strat_by_CDS[mutation_per_day_strat_by_CDS["cds"]=="S"]
# plt.plot(to_plot["date"], to_plot["cumsum"])
# locs, labels = plt.xticks()
# plt.xticks(locs, labels, rotation=45, ha='right')


if __name__ == "__main__":
	main()