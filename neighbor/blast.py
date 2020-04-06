import subprocess
in_fn = "../data/gisaid/EPI_ISL_410218.fasta"

def blast(in_fn):

	op_sys = subprocess.run("uname", capture_output=True).stdout
	if op_sys == b"Linux\n":
		makeblastdb = ""
	elif op_sys == b"Darwin\n":
		blastn = "ext/ncbi-blast-2.10.0+-x64-macosx/bin/blastn"
		db_prefix = "ext/ncbi-blast-2.10.0+-x64-macosx/blastdb/gisaid"
	else:
		raise Exception("Coviz only supports Linux and MacOS!")

	cmd = f"{blastn} -db {db_prefix} -query {in_fn} -outfmt 6 -num_alignments 5".split()
	output = subprocess.run(cmd, capture_output=True)
	query_id, subj_id = output.stdout.decode("utf-8").strip().split("\n")[0].split("\t")[:2]

	return subj_id

# Test this code tomorrow. 