import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def get_mutation_count(in_fn):
	with open(in_fn, "r") as f:
		mut_count = 0
		for line in f:
			if line.startswith("#"):
				continue
			mut_count += 1
	return mut_count


def get_mutation_count_from_dir(in_dir):
	sample_list = []
	count_list = []
	for vcf_fn in glob.glob(f"{in_dir}/*.vcf"):
		print(vcf_fn)
		sample = os.path.basename(vcf_fn).replace(".vcf", "")
		sample_list.append(sample)
		count = get_mutation_count(vcf_fn)
		count_list.append(count)
	assert len(sample_list) == len(count_list)
	df = pd.DataFrame({
		"Sample": sample_list,
		"Mutation Count": count_list
		})
	return df


def plot_mut_per_sample(count, out_fn, vline=200):
	fig, ax = plt.subplots(1,1)
	ax = sns.distplot(count, kde=False, rug=False)
	ax.set_yscale("log")
	ax.set_xlabel("Mutation Count")
	ax.set_ylabel("Frequency")
	plt.axvline(x=vline, color="red", linestyle="--")
	plt.text(x=vline+2, y=100, s=f"count = {vline}", color="red")
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	fig.set_size_inches(4, 3)
	fig.tight_layout()
	fig.savefig(out_fn)


def main():
	in_dir = "../processed_data/fasta2vcf/fasta2vcf/"
	out_dir = "../processed_data/fasta2vcf/mut_per_sample/"
	os.makedirs(out_dir, exist_ok=True)
	df = get_mutation_count_from_dir(in_dir)
	count = df["Mutation Count"].tolist()
	out_fn = f"{out_dir}/mut_per_sample.pdf"
	plot_mut_per_sample(count, out_fn, vline=100)

	count_2 = [x for x in count if x < 200]
	out_fn = f"{out_dir}/mut_per_sample_lt_200.pdf"
	plot_mut_per_sample(count_2, out_fn, vline=50)

if __name__ == "__main__":
	main()
