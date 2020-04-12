import click
import os
import glob
import gzip
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def get_mutation_count(in_fn):
	mut_count = 0

	if ".gz" in in_fn:
		with gzip.open(in_fn, "r") as f:
			for line in f:
				line = line.decode("utf-8")
				if line.startswith("#"):
					continue
				mut_count += 1
	else:
		with open(in_fn, "r") as f:
			for line in f:
				if line.startswith("#"):
					continue
				mut_count += 1

	return mut_count



def get_mutation_count_from_dir(in_dir, out_fn):
	sample_list = []
	count_list = []
	for vcf_fn in glob.glob(f"{in_dir}/*.vcf.gz"):
		print(vcf_fn)
		sample = os.path.basename(vcf_fn).replace(".vcf.gz", "")
		sample_list.append(sample)
		count = get_mutation_count(vcf_fn)
		count_list.append(count)
	assert len(sample_list) == len(count_list)
	mutation_count = pd.DataFrame({
		"Sample": sample_list,
		"Mutation Count": count_list
		})
	mutation_count.to_csv(out_fn, index=False, sep="\t")
	return mutation_count


def write_filtered_vcf_list(in_dir, out_fn, mutation_count, cutoff=150):
	count = 0
	with open(out_fn, "w") as fout:
		for index, row in mutation_count.iterrows():
			sample = row["Sample"]
			mutation_count = row["Mutation Count"]
			if mutation_count <= cutoff:
				fout.write(f"{in_dir}/{sample}.vcf.gz\n")
				count += 1
	print(f"{count} samples passed cutoff = {cutoff}.")


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


@click.command()
@click.option("--in_dir", "-i", type=str, help="VCF directory.")
@click.option("--out_dir", "-o", type=str, help="Output directory.")
@click.option("--cutoff", "-c", type=int, help="Cutoff to filter VCF files. VCF files with more records than cutoff will be removed.")
def main(in_dir, out_dir, cutoff):
	os.makedirs(out_dir, exist_ok=True)
	mutation_count = get_mutation_count_from_dir(in_dir, f"{out_dir}/mutation_count.tsv")
	write_filtered_vcf_list(in_dir, f"{out_dir}/filtered_vcf.txt", mutation_count, cutoff)
	# count = df["Mutation Count"].tolist()
	# out_fn = f"{out_dir}/mut_per_sample.pdf"
	# plot_mut_per_sample(count, out_fn, vline=100)

	# count_2 = [x for x in count if x < 200]
	# out_fn = f"{out_dir}/mut_per_sample_lt_200.pdf"
	# plot_mut_per_sample(count_2, out_fn, vline=50)

if __name__ == "__main__":
	main()
