import pandas as pd
from Bio import SeqIO
import os

fasta_fn = "../processed_data/phylogenetic/filter_distant_seq/keep.fasta"
meta_fn = "../data/gisaid/metadata/metadata.tsv"
out_dir = "../processed_data/phylogenetic/sample_seq/"
os.makedirs(out_dir, exist_ok=True)

def read_fasta_headers(fasta_fn):
	fasta_header_list = []
	for record in SeqIO.parse(fasta_fn, "fasta"):
		fasta_header_list.append(record.id)
	return fasta_header_list


def read_meta(meta_fn):
	meta = pd.read_table(meta_fn)
	return meta


def subset_meta_to_fasta_header(meta, fasta_header_list):
	meta = meta[meta["strain"].isin(fasta_header_list)]
	return meta


def get_region_count(meta):
	region_count = meta.groupby("region").agg(region_count=("region", "count"))
	return region_count


def sample_meta(meta, by, n):
	sampled = meta.groupby(by, group_keys=False).apply(lambda x: x.sample(n))
	return sampled


def filter_fasta(fasta_fn, keep_list, out_fn):
	with open(out_fn, "w") as f:
		n = 0
		for record in SeqIO.parse(fasta_fn, "fasta"):
			if record.id in keep_list:
				SeqIO.write(record, f, "fasta")
				n += 1
				print(f"{n} records written.")


def main():
	fasta_header_list = read_fasta_headers(fasta_fn)
	meta = read_meta(meta_fn)
	meta = subset_meta_to_fasta_header(meta, fasta_header_list)
	region_count = get_region_count(meta)
	sample = sample_meta(meta, by="region", n=30)
	keep_list = sample["strain"].tolist()
	out_fn = f"{out_dir}/sample.fasta"
	filter_fasta(fasta_fn, keep_list, out_fn)


if __name__ == "__main__":
	main()