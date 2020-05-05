import os
import subprocess
import glob
from Bio import SeqIO
import click
import sys
sys.path.append(".")
from utils import DefaultOrderedDict


def get_IDs_from_align_file(align_fn):
	ID = []
	for record in SeqIO.parse(align_fn, "fasta"):
		ID.append(record.id)
	if len(ID) > 2:
		print("Alignment file have >2 sequences. Using the first 2 sequences.")
	return ID[0], ID[1]


def align(ref_fasta, qry_fasta, out_prefix):
	mafft_in_fn = f"{out_prefix}.fasta"
	with open(mafft_in_fn, "w") as fout:
		for fn in [ref_fasta, qry_fasta]:
			with open(fn, "r") as fin:
				for line in fin:
					fout.write(line.strip() + "\n")

	op_sys = subprocess.check_output("uname")
	if op_sys == b"Linux\n":
		mafft = "ext/mafft-linux64/mafft.bat"
	elif op_sys == b"Darwin\n":
		mafft = "ext/mafft-mac/mafft.bat"
	else:
		raise Exception("Coviz only supports Linux and MacOS!")

	mafft_out_fn = f"{out_prefix}.ali"
	cmd = f"{mafft} {mafft_in_fn} > {mafft_out_fn}"
	output = subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL)

	return mafft_out_fn


def parse_align(align_fn):
	ref = ""
	qry = ""
	n_record = 0
	with open(align_fn, "r") as f:
		for line in f:
			if line.startswith(">"):
				if n_record == 0:
					is_ref = True
					n_record += 1
				else:
					is_ref = False
			else:
				if is_ref:
					ref += line.strip()
				else:
					qry += line.strip()

	return ref, qry


def align2variant(ref, qry):
	assert len(ref) == len(qry)
	ref_coord = 0
	qry_coord = 0
	ref_variant = DefaultOrderedDict(str)
	qry_variant = DefaultOrderedDict(str)
	r0 = ""
	q0 = ""

	for r, q in zip(ref,qry):
		r = r.upper().replace("U", "T")
		q = q.upper().replace("U", "T")

		if r != "-":
			ref_coord += 1
		if q != "-":
			qry_coord = ref_coord
		if r != "-" and q != "-":
			r0 = r
			q0 = q

		if r == q:
			pass
		elif r == "n" or q == "n":
			pass
		elif r == "-":
			if ref_variant[ref_coord] == "":
				ref_variant[ref_coord] = r0
				qry_variant[ref_coord] = q0
			ref_variant[ref_coord] += r
			qry_variant[ref_coord] += q
		elif q == "-":
			if ref_variant[qry_coord] == "":
				ref_variant[qry_coord] = r0
				qry_variant[qry_coord] = q0
			ref_variant[qry_coord] += r
			qry_variant[qry_coord] += q
		elif r != q:
			ref_variant[ref_coord] = r
			qry_variant[ref_coord] = q
		else:
			raise Exception("Error! Please email your sequence to jollier.liu@gmail.com")
	return ref_variant, qry_variant


def save_vcf(ref_variant, qry_variant, qry_name, out_fn):
	assert len(ref_variant) == len(qry_variant)
	with open(out_fn, "w") as f:
		f.write('##fileformat=VCFv4.2\n')
		f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
		f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
		f.write('##contig=<ID=NC_045512.2,length=29903,assembly=NC_045512.2>\n')
		f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{qry_name}\n")
		for coord in ref_variant.keys():
			if coord == 0: # skip coord = 0
				continue

			rv = ref_variant[coord].replace("-","")
			qv = qry_variant[coord].replace("-","")

			rv_is_canonical = all([x in ["A","T","G","C"] for x in rv])
			qv_is_canonical = all([x in ["A","T","G","C"] for x in qv])

			if rv_is_canonical and qv_is_canonical:
				f.write(f"NC_045512.2\t{coord}\t.\t{rv}\t{qv}\t.\tPASS\t.\tGT\t1\n")


def filter_polya(vcf_fn):
	print("Removing variants in the Poly-A tail.")
	with open(vcf_fn, "r") as fin, open(f"{vcf_fn}.tmp", "w") as fout:
		for line in fin:
			if line.startswith("#"):
				fout.write(line)
			else:
				split_line = line.strip().split("\t")
				pos = split_line[1]
				ref = split_line[3]
				if pos == "29870" and ref.endswith("AAAAAAAAAA"):
					pass
				else:
					fout.write(line)
	os.remove(vcf_fn)
	os.rename(f"{vcf_fn}.tmp", vcf_fn)


def postprocess_vcf(vcf_fn, ref_fn, compress_vcf):
	'''Normalize, rename ID, bgzip and tabix'''
	print("Normalizing VCF.")
	cmd = f"bcftools norm -f {ref_fn} {vcf_fn} -Ou | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -Ov -o {vcf_fn}"
	output1 = subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL)
	filter_polya(vcf_fn)

	if compress_vcf:
		print("Compressing and indexing VCF.")
		cmd = f"bgzip {vcf_fn}"
		output2 = subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL)
		cmd = f"tabix -p vcf {vcf_fn}.gz"
		output3 = subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL)

def cleanup(out_prefix, compress_vcf):
	os.remove(f"{out_prefix}.ali")
	os.remove(f"{out_prefix}.fasta")


def align2vcf(align_fn, ref_fn, qry_id, out_prefix, compress_vcf, clean_up, verbose=False):
	if verbose: print("Parsing alignment...")
	ref, qry = parse_align(align_fn)

	if verbose: print("Making VCF...")
	ref_variant, qry_variant = align2variant(ref, qry)

	out_fn = f"{out_prefix}.vcf"
	if verbose: print(f"Saving VCF to {out_fn}")
	save_vcf(ref_variant, qry_variant, qry_id, out_fn)

	if verbose: print(f"Normalize, update ID, and index.")
	postprocess_vcf(out_fn, ref_fn, compress_vcf)

	if clean_up:
		if verbose: print(f"Removing .ali, .fasta, and .vcf files.")
		cleanup(out_prefix, compress_vcf)


def fasta2vcf(fasta_fn, ref_fn, align_fn, out_dir, compress_vcf, clean_up, verbose):
	os.makedirs(out_dir, exist_ok=True)

	# If user set align_fn option.
	if align_fn:
		ref_id, qry_id = get_IDs_from_align_file(align_fn)
		out_prefix = f"{out_dir}/{qry_id}"

		align2vcf(align_fn, ref_fn, qry_id, \
			out_prefix, compress_vcf, clean_up, verbose)

	# If align_fn is empty, start from fasta_fn + ref_fn.
	else:
		for qry_record in SeqIO.parse(fasta_fn, "fasta"):
			qry_id = qry_record.id
			print(f"Query Name: {qry_id}")

			out_prefix = f"{out_dir}/{qry_id}"
			if os.path.exists(f"{out_prefix}.vcf.gz"):
				print(f"{out_prefix}.vcf.gz already exists!")
				continue

			if verbose: print("Aligning...")
			with open("qry.fasta", "w") as f:
				SeqIO.write(qry_record, f, "fasta")
			align_fn = align(ref_fn, "qry.fasta", out_prefix)
			os.remove("qry.fasta")

			align2vcf(align_fn, ref_fn, qry_id, \
				out_prefix, compress_vcf, clean_up, verbose)

@click.command()
@click.option("--fasta_fn", "-f", type=str, help="Input FASTA file.")
@click.option("--ref_fn", "-r", type=str, help="Reference FASTA file.")
@click.option("--align_fn", "-a", type=str, help="Alignment FASTA file. If align_fn is set, program will prioritize align_fn over fasta_fn + ref_fn. Note: align_fn must contain only 1 reference followed by only 1 query sequence.")
@click.option("--out_dir", "-o", type=str, help="Output directory.")
@click.option("--compress_vcf", type=bool, default=True, help="Whether to compress VCF file. Default to True.")
@click.option("--clean_up", type=bool, default=True, help="Whether to clean up files such as .ali and .fasta.")
@click.option("--verbose", "-v", is_flag=True, default=False)
def main(fasta_fn, ref_fn, align_fn, out_dir, compress_vcf, clean_up, verbose):
	print("################")
	print("# FASTA to VCF #")
	print("################")
	print(f"FASTA: {fasta_fn}")
	print(f"Reference: {ref_fn}")
	print(f"Align: {align_fn}")
	print(f"Output: {out_dir}")
	print(f"Compress VCF: {compress_vcf}")
	print(f"Clean up: {clean_up}")
	print(f"Verbose: {verbose}")

	fasta2vcf(fasta_fn, ref_fn, align_fn, out_dir, compress_vcf, clean_up, verbose)

if __name__ == "__main__":
	main()