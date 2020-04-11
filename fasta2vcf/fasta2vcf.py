import os
import subprocess
import glob
from collections import defaultdict
from collections import OrderedDict, Callable

class DefaultOrderedDict(OrderedDict):
	# Source: http://stackoverflow.com/a/6190500/562769
	def __init__(self, default_factory=None, *a, **kw):
		if (default_factory is not None and
		   not isinstance(default_factory, Callable)):
			raise TypeError('first argument must be callable')
		OrderedDict.__init__(self, *a, **kw)
		self.default_factory = default_factory

	def __getitem__(self, key):
		try:
			return OrderedDict.__getitem__(self, key)
		except KeyError:
			return self.__missing__(key)

	def __missing__(self, key):
		if self.default_factory is None:
			raise KeyError(key)
		self[key] = value = self.default_factory()
		return value

	def __reduce__(self):
		if self.default_factory is None:
			args = tuple()
		else:
			args = self.default_factory,
		return type(self), args, None, None, self.items()

	def copy(self):
		return self.__copy__()

	def __copy__(self):
		return type(self)(self.default_factory, self)

	def __deepcopy__(self, memo):
		import copy
		return type(self)(self.default_factory,
						  copy.deepcopy(self.items()))

	def __repr__(self):
		return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
											   OrderedDict.__repr__(self))


def get_sequence_length(fasta):
	seq_length = 0
	with open(fasta) as f:
		for line in f:
			if line.startswith(">"):
				continue
			seq_length += len(line.strip())
	return seq_length

def align(ref_fasta, qry_fasta, out_prefix):
	print("Aligning reference and query sequences...")

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
	print("Parsing alignment...")
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


def main():
	in_dir = "../data/gisaid/"
	ref_fasta = "../data/reference/NC_045512.2.fasta"
	out_dir = "../processed_data/fasta2vcf/fasta2vcf"
	os.makedirs(out_dir, exist_ok=True)
	print(os.getcwd())
	for qry_fasta in glob.glob(f"{in_dir}/*.fasta"):
		print(qry_fasta)

		seq_length = get_sequence_length(qry_fasta)
		if seq_length < 29000:
			print("Too short. Skip!")
			continue
		
		qry_name = os.path.basename(qry_fasta).replace('.fasta', "")
		print(f"Query Name: {qry_name}")

		out_prefix = f"{out_dir}/{qry_name}"

		print("Aligning...")
		align_fn = align(ref_fasta, qry_fasta, out_prefix)
		ref, qry = parse_align(align_fn)

		print("Making VCF...")
		ref_variant, qry_variant = align2vcf(ref, qry)

		out_fn = f"{out_prefix}.vcf"
		print(f"Saving VCF to {out_fn}")
		save_vcf(ref_variant, qry_variant, qry_name, out_fn)


if __name__ == "__main__":
	main()