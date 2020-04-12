import os
import subprocess
import glob
from Bio import SeqIO
import click
import sys
sys.path.append(".")
from utils import DefaultOrderedDict

def get_sequence_length(fasta):
	seq_length = 0
	with open(fasta) as f:
		for line in f:
			if line.startswith(">"):
				continue
			seq_length += len(line.strip())
	return seq_length

def align(ref_fasta, qry_fasta, out_prefix):
	mafft_in_fn = f"{out_prefix}.fasta"
	with open(mafft_in_fn, "w") as fout:
		for fn in [ref_fasta, qry_fasta]:
			with open(fn, "r") as fin:
				for line in fin:
					fout.write(line.strip() + "\n")


	mafft_out_fn = f"{out_prefix}.ali"
	cmd = f"mafft {mafft_in_fn} > {mafft_out_fn}"
	output = subprocess.run(cmd, shell=True, capture_output=True)

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


def align2vcf(ref, qry):
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
		f.write('##contig=<ID=1,length=29903,assembly=NC_045512.2>\n')
		f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{qry_name}\n")
		for coord in ref_variant.keys():
			if coord == 0: # skip coord = 0
				continue

			rv = ref_variant[coord].replace("-","")
			qv = qry_variant[coord].replace("-","")

			rv_is_canonical = all([x in ["A","T","G","C"] for x in rv])
			qv_is_canonical = all([x in ["A","T","G","C"] for x in qv])

			if rv_is_canonical and qv_is_canonical:
				id_ = f"1_{coord}_{rv[:5]}_{qv[:5]}"
				f.write(f"1\t{coord}\t{id_}\t{rv}\t{qv}\t.\tPASS\t.\tGT\t1\n")


def postprocess_vcf(vcf_fn, ref_fn):
	'''Normalize, rename ID, bgzip and tabix'''
	cmd = f"bcftools norm -f {ref_fn} {vcf_fn} -Ou | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -Oz -o {vcf_fn}.gz"
	output1 = subprocess.run(cmd, shell=True, capture_output=True)
	cmd = f"tabix -p vcf {vcf_fn}.gz"
	output2 = subprocess.run(cmd, shell=True, capture_output=True)


def cleanup(out_prefix):
	os.remove(f"{out_prefix}.ali")
	os.remove(f"{out_prefix}.fasta")
	os.remove(f"{out_prefix}.vcf")


@click.command()
@click.option("--fasta_fn", "-f", type=str, help="Input FASTA file.")
@click.option("--ref_fn", "-r", type=str, help="Reference FASTA file.")
@click.option("--out_dir", "-o", type=str, help="Output directory.")
@click.option("--verbose", "-v", is_flag=True, default=False)
def main(fasta_fn, ref_fn, out_dir, verbose):
	print("################")
	print("# FASTA to VCF #")
	print("################")
	print(f"FASTA: {fasta_fn}")
	print(f"Reference: {ref_fn}")
	print(f"Output: {out_dir}")

	os.makedirs(out_dir, exist_ok=True)

	for qry_record in SeqIO.parse(fasta_fn, "fasta"):
		qry_id = qry_record.id
		print(f"Query Name: {qry_id}")

		out_prefix = f"{out_dir}/{qry_id}"

		if verbose: print("Aligning...")
		with open("qry.fasta", "w") as f:
			SeqIO.write(qry_record, f, "fasta")
		align_fn = align(ref_fn, "qry.fasta", out_prefix)
		os.remove("qry.fasta")

		if verbose: print("Parsing alignment...")
		ref, qry = parse_align(align_fn)


		if verbose: print("Making VCF...")
		ref_variant, qry_variant = align2vcf(ref, qry)

		out_fn = f"{out_prefix}.vcf"
		if verbose: print(f"Saving VCF to {out_fn}")
		save_vcf(ref_variant, qry_variant, qry_id, out_fn)

		if verbose: print(f"Normalize, update ID, and index.")
		postprocess_vcf(out_fn, ref_fn)

		if verbose: print(f"Removing .ali, .fasta, and .vcf files.")
		cleanup(out_prefix)

if __name__ == "__main__":
	main()