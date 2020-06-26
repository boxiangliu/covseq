from Bio import SeqIO
import time
import subprocess
from multiprocessing import Pool
import os

ref_fasta_fn = "data/NC_045512.2.fasta"
qry_fasta_fn = "../data/gisaid/fasta/sequences.fasta"
out_dir = "../processed_data/phylogenetic/filter_distant_seq/"
os.makedirs(out_dir, exist_ok=True)


def write_fasta(ref_fasta_fn, qry_fasta_fn, out_dir):
	out_fn_list = []
	ref_fasta = SeqIO.read(ref_fasta_fn, "fasta")
	for i, qry_fasta in enumerate(SeqIO.parse(qry_fasta_fn, "fasta")):
		out_fn = f"{out_dir}/{i:05}.fasta"
		out_fn_list.append(out_fn)
		with open(out_fn, "w") as f:
			SeqIO.write(ref_fasta, f, "fasta")
			SeqIO.write(qry_fasta, f, "fasta")
	return out_fn_list


def pairwise_alignment(fasta_fn, align_fn=None, verbose=False):
	if align_fn is None:
		align_fn = fasta_fn.replace(".fasta", ".ali")
	cmd = f"mafft --retree 1 --maxiterate 0 {fasta_fn} > {align_fn}"
	if verbose:
		print(f"Command: {cmd}")
	output = subprocess.run(cmd, shell=True, capture_output=True)
	return align_fn


def get_hamming_distance(fasta_fn):
	records = SeqIO.parse(fasta_fn, "fasta")
	ref_fasta = next(records)
	qry_fastt = next(records)
	hamming = sum([x != y for x, y in zip(str(ref_fasta.seq), str(qry_fastt.seq))])
	return hamming


def filter_seq_by_hamming_dist(qry_fasta_fn, hamming, out_dir):
	records = SeqIO.parse(qry_fasta_fn, "fasta")
	keep_fn = f"{out_dir}/keep.fasta"
	remove_fn = f"{out_dir}/remove.fasta"
	hamming_fn = f"{out_dir}/hamming.tsv"

	with open(keep_fn, "w") as f_keep, \
		open(remove_fn, "w") as f_remove, \
		open(hamming_fn, "w") as f_hamming:
		for i, (h, qry_fasta) in enumerate(zip(hamming, records)):
			if h < 300:
				SeqIO.write(qry_fasta, f_keep, "fasta")
			else:
				SeqIO.write(qry_fasta, f_remove, "fasta")
			f_hamming.write(f"{i}\t{qry_fasta.id}\t{h}\n")


def main():
	print("Start filter_distant_seq.py")
	start = time.time()

	print("Write FASTA for pairwise alignment.")
	fasta_fn_list = write_fasta(ref_fasta_fn, qry_fasta_fn, out_dir)

	write_time = time.time()
	write_duration = write_time - start
	print(f"Writing took {write_duration} seconds.")


	print("Pairwise alignment.")
	with Pool(40) as p:
		align_fn_list = p.map(pairwise_alignment, fasta_fn_list)

	align_time = time.time()
	align_duration = align_time - write_time
	print(f"Pairwise alignment took {align_duration} seconds.")


	print("Hamming distance.")
	with Pool(40) as p:
		hamming = p.map(get_hamming_distance, align_fn_list)

	hamming_time = time.time()
	hamming_duration = hamming_time - align_time
	print(f"Hamming distance calculation took {hamming_duration} seconds.")


	print("Filter by Hamming distance.")
	filter_seq_by_hamming_dist(qry_fasta_fn, hamming, out_dir)

	filter_time = time.time()
	filter_duration = filter_time - hamming_time
	print(f"Filtering took {filter_duration} seconds.")

	finish_time = time.time()
	finish_duration = finish_time - start
	print(f"Total time: {finish_duration} seconds.")


# Writing took 60.859065771102905 seconds.
# Pairwise alignment took 592.5926842689514 seconds.
# Hamming distance calculation took 9.643347263336182 seconds.
# Filtering took 12.368592262268066 seconds. 
# Total time: 671.4485294818878.
