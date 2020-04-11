'''
Filter out samples with too many mutations.
'''
import os
import glob
import shutil

def get_mutation_count(in_fn):
	with open(in_fn, "r") as f:
		mut_count = 0
		for line in f:
			if line.startswith("#"):
				continue
			mut_count += 1
	return mut_count



def main():
	in_dir = "../processed_data/fasta2vcf/fasta2vcf/"
	out_dir = "../processed_data/fasta2vcf/filter_samples/"
	threshold = 50
	os.makedirs(out_dir, exist_ok=True)

	for vcf_fn in glob.glob(f"{in_dir}/*.vcf"):
		print(vcf_fn)
		name = os.path.basename(vcf_fn)
		count = get_mutation_count(vcf_fn)
		print(f"Number of records: {count}")
		if count >= threshold:
			continue
		shutil.copyfile(f"{vcf_fn}", f"{out_dir}/{name}")


if __name__ == "__main__":
	main()