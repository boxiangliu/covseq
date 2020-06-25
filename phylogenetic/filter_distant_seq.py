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


def pairwise_alignment(fasta_fn, align_fn=None):
	if align_fn is None:
		align_fn = fasta_fn.replace(".fasta", ".ali")
	cmd = f"mafft --retree 1 --maxiterate 0 {fasta_fn} > {align_fn}"
	print(f"Command: {cmd}")
	subprocess.run(cmd, shell=True)
	return align_fn


def get_hamming_distance(fasta_fn):
	records = SeqIO.parse(fasta_fn, "fasta")
	ref_fasta = next(records)
	qry_fastt = next(records)
	hamming = sum([x != y for x, y in zip(str(ref_fasta.seq), str(qry_fastt.seq))])
	return hamming


def filter_seq_by_hamming_dist(qry_fasta_fn, hamming, out_fn):
	records = SeqIO.parse(qry_fasta_fn, "fasta")
	with open(out_fn, "w") as f:
		for h, qry_fasta in zip(hamming, records):
			if h < 300:
				SeqIO.write(qry_fasta, f, "fasta")


start = time.time()

fasta_fn_list = write_fasta(ref_fasta_fn, qry_fasta_fn, out_dir)

write_time = time.time()
write_duration = write_time - start
print(f"Writing took {write_duration} seconds.")

with Pool(40) as p:
	align_fn_list = p.map(pairwise_alignment, fasta_fn_list)

align_time = time.time()
align_duration = align_time - write_time
print(f"Pairwise alignment took {align_duration} seconds.")


with Pool(40) as p:
	hamming = p.map(get_hamming_distance, align_fn_list)

hamming_time = time.time()
hamming_duration = hamming_time - align_time
print(f"Hamming distance calculation took {hamming_duration} seconds.")


out_fn = f"{out_dir}/filtered.fasta"
filter_seq_by_hamming_dist(qry_fasta_fn, hamming, out_fn)

filter_time = time.time()
filter_duration = filter_time - hamming_time
print(f"Filtering took {filter_duration}.")

finish_time = time.time()
finish_duration = finish_time - start
print(f"Total time: {finish_duration}.")
