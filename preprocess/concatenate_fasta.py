import glob
import os
import gzip
from Bio import SeqIO
in_dir = "../data/"
sub_dirs = ["gisaid", "ncbi", "cngb", "embl"]
out_fn = "../data/aggregated/fasta/raw.fasta"
cgnb_metadata = "../data/cngb/metadata/CNGBdb_VirusDIP_excel20200411_all(24)_57b4ce53c4d6c49c978596677a112211.csv"

print("#########################")
print("# Concatenate all FASTA #")
print("#########################")
print(f"Input: {in_dir}")
print(f"Output: {out_fn}")


def get_cngb_id_map(in_fn):
	id_map = {}
	with open(in_fn) as f:
		for line in f:
			if "Sequence ID" in line:
				continue
			split_line = line.split(",")
			fasta_fn = os.path.basename(split_line[-2]).replace(".gz", "")
			ID = split_line[0]
			id_map[fasta_fn] = ID
	return id_map


def concatenate_fasta(in_dir, sub_dirs, out_fn, cngb_id_map):
	counter = 0
	with open(out_fn, "w") as fout:
		for sd in sub_dirs:
			wd = f"{in_dir}/{sd}/fasta/"
			print(f"Working on directory {wd}")
			for fn in glob.glob(f"{wd}/*.fasta"):
				for record in SeqIO.parse(fn, "fasta"):
					if sd == "gisaid":
						seq_header = record.description.split("|")[1].strip()
					elif sd == "ncbi":
						seq_header = record.description.split("|")[0].strip()
					elif sd == "embl":
						seq_header = record.description.split("|")[1].strip()
					elif sd == "cngb":
						bn = os.path.basename(fn)
						seq_header = cngb_id_map[bn]
					record.id = record.name = record.description = seq_header
					SeqIO.write(record, fout, "fasta")
					counter += 1
	print(f"Combined {counter} sequences!")


def main():
	cngb_id_map = get_cngb_id_map(cgnb_metadata)
	concatenate_fasta(in_dir, sub_dirs, out_fn, cngb_id_map)

if __name__ == "__main__":
	main()
