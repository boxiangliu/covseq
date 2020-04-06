import subprocess
import glob

in_dir = "/Users/boxiang/Documents/work/Baidu/projects/viraviz/data/gisaid/"
out_fn = "data/gisaid.fasta"
def prepare_input_file(in_dir, out_fn):
	print("Making input fasta file...")
	with open(out_fn, "w") as fout:
		for fn in glob.glob(f"{in_dir}/*.fasta"):
			print(fn)
			with open(fn, "r") as fin:
				for line in fin:
					line = line.strip()
					if line.startswith(">"):
						header = line.split("|")[1]
						fout.write(f">{header}\n")
					else:
						fout.write(f"{line}\n")



def make_blast_db(in_fn):
	print("Making blast database...")
	op_sys = subprocess.run("uname", capture_output=True).stdout
	if op_sys == b"Linux\n":
		makeblastdb = ""
	elif op_sys == b"Darwin\n":
		makeblastdb = "ext/ncbi-blast-2.10.0+-x64-macosx/bin/makeblastdb"
		db_dir = "ext/ncbi-blast-2.10.0+-x64-macosx/blastdb/"

	else:
		raise Exception("Coviz only supports Linux and MacOS!")

	cmd = f"{makeblastdb} -in {in_fn} -dbtype nucl -parse_seqids -out {db_dir}/gisaid -title gisaid"
	output = subprocess.run(cmd, shell=True, capture_output=True)

	return output

if __name__ == "__main__":
	prepare_input_file(in_dir, out_fn)
	output = make_blast_db(out_fn)