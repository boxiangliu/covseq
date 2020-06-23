import os
from pysam import VariantFile
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

vcf_fn = "../data/aggregated/vcf/merged/filtered.vcf.gz"
out_dir = "../processed_data/mutation_distribution/sample_per_mutation/"
os.makedirs(out_dir, exist_ok=True)

def tally_mutations_per_sample(vcf):
	record = next(vcf)
	total_samples = len(record.samples.items())
	vcf.reset()
	mut_count_dict = defaultdict(int)
	for i, record in enumerate(vcf):
		print(f"Record: {i}")
		mut_count = 0
		ref = record.ref
		for samples in record.samples.items():
			alleles = samples[1].alleles[0]
			if alleles != ref:
				mut_count += 1
		mut_count_dict[record.pos] += mut_count

	pos_list = []
	mut_count_list = []
	for pos, count in mut_count_dict.items():
		pos_list.append(pos)
		mut_count_list.append(count)

	proportion = [x/total_samples for x in mut_count_list]
	df = pd.DataFrame({"pos":pos_list, "count":mut_count_list, "proportion": proportion})
	df.sort_values(by="pos", inplace=True)

	return df


def plot_mutations_per_sample(sample_per_mut, out_fn):
	fig, ax = plt.subplots(1,1)
	ax = plt.bar(x=sample_per_mut["pos"], height=sample_per_mut["proportion"])
	plt.xticks(ticks=range(0,30001,2000), labels=range(0,30001,2000))
	font = {'family': 'sans-serif', 'color':  'black', 'weight': 'bold', 'size': 12}
	plt.xlabel("Genome coordinate", fontdict=font)
	plt.ylabel("Proportion mutated", fontdict=font)
	fig.set_size_inches(12, 3)
	fig.tight_layout()
	fig.savefig(out_fn)


def main():
	vcf = VariantFile(vcf_fn)
	sample_per_mut = tally_mutations_per_sample(vcf)
	sample_per_mut.to_csv(f"{out_dir}/sample_per_mut.tsv", sep="\t", index=False)
	plot_mutations_per_sample(sample_per_mut, f"{out_dir}/sample_per_mut.pdf")


if __name__ == "__main__":
	main()