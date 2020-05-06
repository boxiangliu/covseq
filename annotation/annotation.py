import click
import os
import subprocess
from collections import defaultdict, OrderedDict
import numpy as np
import pandas as pd
from Bio import SeqIO, SeqFeature, Seq
import sys
sys.path.append(".")
from vcf.fasta2vcf import fasta2vcf
from utils import DefaultOrderedDict, VCF
from snpEff.parse_snpEff import HEADERS, parse_snpEff
from werkzeug.utils import secure_filename

class Error():
	def __init__(self, id_, msg, qryid):
		self.id = id_
		self.msg = msg
		self.qryid = qryid


def read_fasta(fasta_fn):
	fasta = list(SeqIO.parse(fasta_fn, "fasta"))
	for i, x in enumerate(fasta):
		x.id = secure_filename(x.id)

		print("Replacing U's with T's.")
		x.seq = str(x.seq).upper().replace("U","T")

		non_ATGC = 0
		for nuc in x.seq:
			if nuc not in "ATGC":
				non_ATGC += 1
		
		if non_ATGC/len(x.seq) > 0.05:
			fasta[i] = Error(0, "ERROR: sequence contains more than 5% non-A[T/U]GC letters.", x.id)

	return fasta


def read_ref_genbank(gbk_fn):
	return next(SeqIO.parse(gbk_fn, "genbank"))


class Annotation():
	def __init__(self, qry, ref_gbk, out_dir, verbose):
		'''
		qry: a SeqIO FASTA object.
		ref_gbk: a SeqIO genbank object.
		out_dir: output directory.
		'''
		self.qry = qry
		self.ref_gbk = ref_gbk
		self.out_dir = out_dir
		self.verbose = verbose

	def run(self):
		print("Start annotating the query sequence.")
		print(f"Query ID: {self.qry.id}")
		print(f"Ref. ID: {self.ref_gbk.id}")
		if len(self.qry.seq) < 25000:
			print("Query sequence should have at least 25,000 nucleotides! Aborting.")
			return 

		self.align(self.ref_gbk, self.qry, self.out_dir)
		self.parse_align(self.align_fn)
		# self.get_mutation(self.nt_df)
		# self.transfer_feature(self.ref_gbk, self.nt_df, self.verbose)
		self.get_orf(self.ref_gbk, self.nt_df)


	def align(self, ref, qry, out_dir):
		print("Aligning reference and query sequences.")
		work_dir = f"{out_dir}/{qry.id}"
		os.makedirs(work_dir, exist_ok=True)
		
		op_sys = subprocess.check_output("uname")
		if op_sys == b"Linux\n":
			mafft = "ext/mafft-linux64/mafft.bat"
		elif op_sys == b"Darwin\n":
			mafft = "ext/mafft-mac/mafft.bat"
		else:
			raise Exception("Coviz only supports Linux and MacOS!")

		mafft_in_fn = f"{work_dir}/{qry.id}.fasta"
		with open(mafft_in_fn, "w") as f:
			f.write(f">{ref.id}\n")
			f.write(f"{str(ref.seq)}\n")
			f.write(f">{qry.id}\n")
			f.write(f"{str(qry.seq)}\n")

		mafft_out_fn = f"{work_dir}/{qry.id}.ali"
		cmd = f"{mafft} {mafft_in_fn} > {mafft_out_fn}"
		output = subprocess.check_output(cmd, stderr=subprocess.DEVNULL, shell=True)

		self.align_fn = mafft_out_fn


	def parse_align(self, align_fn):
		print("Parsing alignment.")
		align = defaultdict(list)
		coord = defaultdict(list)
		header = 0
		with open(align_fn, "r") as f:
			for line in f:
				if line.startswith(">"):
					header += 1
					if header == 1:
						strain = "ref"
					else:
						strain = "qry"

					counter = 1
				else:
					for base in line.strip():
						align[strain].append(base)
						if base != "-":
							coord[strain].append(counter)
							counter += 1
						else:
							coord[strain].append(np.nan)

		ref_coord = pd.Series(coord["ref"], dtype="Int64")
		qry_coord = pd.Series(coord["qry"], dtype="Int64")

		nt_df = pd.DataFrame({
			"ref_coord": ref_coord,
			"ref_nt": align["ref"],
			"qry_coord": qry_coord,
			"qry_nt": align["qry"]
			})

		self.nt_df = nt_df 


	def get_mutation(self, nt_df):
		print("Finding nucleotide mutations between ref. and query sequences.")
		ref_nt = nt_df["ref_nt"]
		qry_nt = nt_df["qry_nt"]

		nt_df["nt_mut"] = [f"{i},{j}" if i!=j else "same" \
			for i,j in zip(ref_nt, qry_nt)]

		self.nt_df = nt_df


	@staticmethod
	def init_protein_annotation(location):
		start = location.start
		end = location.end
		protein_anno = OrderedDict()
		for i in range(start, end):
			'''Each entry in protein_anno has 4 elements:
				1) amino acid (default: "", meaning missing value)
				2) frame (default: -1, meaning missing value)
				3) ribosomal slippage (default: 0, meaning no shift)
				4) RNA editing (default: "", meaning no editing)
			'''
			protein_anno[i] = OrderedDict({"aa": "", "frame": -1, "ribo_slip": 0, "rna_edit": ""})
		return protein_anno


	@staticmethod
	def update_protein_annotation(location, seq, protein_anno):

		feature = SeqFeature.SeqFeature(location)
		protein = feature.translate(seq, cds=False)

		start = location.start
		end = location.end

		if (end - start + 1) % 3 == 1:
			end -= 1
		elif (end - start + 1) % 3 == 2:
			end -= 2

		for i in range(start, end+1):

			feature_pos = i - start
			feature_frame = feature_pos % 3
			aa = protein[feature_pos // 3]
			protein_anno[i]["aa"] = aa

			if protein_anno[i]["frame"] != -1:
				frame_diff = protein_anno[i]["frame"] - feature_frame - 3
				protein_anno[i]["ribo_slip"] = frame_diff

			protein_anno[i]["frame"] = feature_frame

		return protein_anno


	@staticmethod
	def get_protein_anno(location, seq):
		protein_anno = Annotation.init_protein_annotation(location)

		if isinstance(location, SeqFeature.FeatureLocation):
			
			protein_anno = Annotation.update_protein_annotation(location, seq, protein_anno)

		elif isinstance(location, SeqFeature.CompoundLocation):
			
			for loc_i in location.__dict__["parts"]:
				protein_anno = Annotation.update_protein_annotation(loc_i, seq, protein_anno)

		else: 

			raise Exception(type(location) + " not defined!")

		return protein_anno


	@staticmethod
	def transfer_simple_location(start, end, fro, to):
		try:
			new_start = int(to[fro == start].values[0])
			new_end = int(to[fro == end].values[0])
			location = SeqFeature.FeatureLocation(new_start, new_end)
		except:
			location = None
		return location


	@staticmethod
	def transfer_location(location, fro, to):

		if isinstance(location, SeqFeature.FeatureLocation):
			# FeatureLocation has only 1 start and 1 end.
			new_location = Annotation.transfer_simple_location(location.start, \
				location.end, fro, to)

		elif isinstance(location, SeqFeature.CompoundLocation):
			# CompoundLocation has multiple locations.
			loc_0 = location.__dict__["parts"][0]
			new_location = Annotation.transfer_simple_location(loc_0.start, \
				loc_0.end, fro, to)

			if new_location is not None:
				for loc_i in location.__dict__["parts"][1:]:
					new_location += Annotation.transfer_simple_location(loc_i.start, \
						loc_i.end, fro, to)

		else:
			raise Exception(type(location) + " not defined!")
		
		return new_location


	@staticmethod
	def dict2df(dict_, colnames, key2col=True):
		
		container = defaultdict(list)
		
		if key2col == True:
		
			for k, v in dict_.items():
				container[colnames[0]].append(k)
				for i, v_i in enumerate(v.values()):
					container[colnames[1+i]].append(v_i)
		
		else:
		
			for k, v in dict_.items():
				for i, v_i in enumerate(v.values()):
					container[colnames[i]].append(v_i)
		
		df = pd.DataFrame(container)
		
		return df


	@staticmethod
	def compare_aa(x, y):
		if x=="" and y=="":
			return ""
		elif x=="":
			return f"-,{y}"
		elif y=="":
			return f"{x},-"
		elif x != y:
			return f"{x},{y}"
		else:
			return "same"


	def get_orf(self, ref_gbk, nt_df, ref_coordinate=False):
		print("Getting ORF annotation.")
		orf = DefaultOrderedDict(list)
		for ref_feature in ref_gbk.features:
			if ref_feature.type == "CDS":
				gene = ref_feature.qualifiers["gene"][0]
				product = ref_feature.qualifiers["product"][0]
				if ref_coordinate:
					location = ref_feature.location
				else:
					location = Annotation.transfer_location(ref_feature.location, \
						nt_df["ref_coord"], nt_df["qry_coord"])
				if location is None:
					continue
				if isinstance(location, SeqFeature.CompoundLocation):
					start = []
					end = []
					for loc_i in location.__dict__["parts"]:
						start.append(str(int(loc_i.start) + 1))
						end.append(str(loc_i.end))
					start = ",".join(start)
					end = ",".join(end)
				else:
					start = str(int(location.start) + 1)
					end = str(location.end)

				rna_length = location.end - location.start + 1
				ribosomal_slippage = "Yes" \
					if "ribosomal_slippage" in ref_feature.qualifiers else "No"
				strand = "+"
				frame = (location.start + 1) % 3
				if frame == 0: frame = 3

				orf["Gene"].append(gene)
				orf["Product"].append(product)
				orf["Start"].append(start)
				orf["End"].append(end)
				orf["Strand"].append(strand)
				orf["Frame"].append(frame)
				orf["RNA_length"].append(rna_length)
				orf["Ribo_Slip"].append(ribosomal_slippage)
		orf = pd.DataFrame(orf)
		self.orf = orf


	def transfer_feature(self, ref_gbk, nt_df, verbose):
		print("Starting transfering CDS fetaures from ref to query.")

		gene_container = OrderedDict()
		for ref_feature in ref_gbk.features:

			if ref_feature.type == "CDS":

				gene = ref_feature.qualifiers["gene"][0].replace(" ","_")
				if verbose: print(f" - ORF: {gene}")
				
				if verbose: print(f" - Annotating reference protein.")
				ref_protein_anno = self.get_protein_anno(ref_feature.location, ref_gbk.seq)

				if verbose: print(f" - Transfering features from reference to query.")
				qry_location = self.transfer_location(ref_feature.location, \
					nt_df["ref_coord"], nt_df["qry_coord"])
				qry_seq = Seq.Seq("".join(nt_df["qry_nt"].tolist()).replace("-",""))
				qry_protein_anno = self.get_protein_anno(qry_location, qry_seq)


				if verbose: print(" - Converting protein annotations into DataFrames.")
				colnames = ["qry_coord", f"{gene}:qry_aa",f"{gene}:frame",
					f"{gene}:ribo_slip",f"{gene}:rna_edit"]
				qry_protein_df = self.dict2df(qry_protein_anno, colnames)
				qry_protein_df["qry_coord"] = qry_protein_df["qry_coord"].astype("Int64")


				colnames = ["ref_coord", f"{gene}:ref_aa", "placeholder1", 
					"placeholder2", "placeholder3"]
				ref_protein_df = self.dict2df(ref_protein_anno, colnames)
				ref_protein_df["ref_coord"] = ref_protein_df["ref_coord"].astype("Int64")
				ref_protein_df.drop(["placeholder1", "placeholder2", "placeholder3"], axis=1, inplace=True)

				gene_container[gene]={"ref": ref_protein_df, "qry": qry_protein_df}


		print("Finding amino acid mutations.")
		merged = nt_df

		for k, v in gene_container.items():
			merged = pd.merge(merged, v["ref"], on="ref_coord", how="outer")
			merged = pd.merge(merged, v["qry"], on="qry_coord", how="outer")
			ref_aa = merged[f"{k}:ref_aa"].fillna("").tolist()
			qry_aa = merged[f"{k}:qry_aa"].fillna("").tolist()

			merged[f"{k}:aa_mut"] = [self.compare_aa(i,j) \
				for i, j in zip(ref_aa, qry_aa)]

		self.anno_df = merged.sort_values("ref_coord")


	def write_sequences(self, out_prefix):
		orf = self.orf
		qry = self.qry

		with open(f"{out_prefix}_nt.fasta", "w") as f_nt, \
			open(f"{out_prefix}_aa.fasta", "w") as f_aa:
			for index, row in orf.iterrows():
				start = row["Start"]
				end = row["End"]
				gene = row["Gene"]

				if "," in start:
					nt_start = int(start.split(",")[0])
					nt_end = int(end.split(",")[-1])
				else:
					nt_start = int(start)
					nt_end = int(end)
				nt_seq = qry.seq[(nt_start-1):nt_end]
				f_nt.write(f">{gene}\n")
				f_nt.write(str(nt_seq) + "\n")

				if "," in start:
					aa_starts = [int(x) for x in start.split(",")]
					aa_ends = [int(x) for x in end.split(",")]
					aa_seq = ""
					for aa_start, aa_end in zip(aa_starts, aa_ends):
						nt_seq = qry.seq[(aa_start-1):aa_end]
						aa_seq += nt_seq.translate()
				else:
					aa_seq = nt_seq.translate()
				f_aa.write(f">{gene}\n")
				f_aa.write(str(aa_seq) + "\n")


def run_snpEff(vcf_fn, out_fn):
	print("Running snpEff.")
	cmd = f"java -jar ext/snpEff/snpEff.jar NC_045512.2 {vcf_fn} > {out_fn}"
	output = subprocess.check_output(cmd, stderr=subprocess.DEVNULL, shell=True)
	os.remove("snpEff_genes.txt")
	os.remove("snpEff_summary.html")


def vcf_intersect_orf(vcf_fn, orf):
	intersect = []
	vcf = VCF(vcf_fn)
	for i, v_row in vcf.rowdata.iterrows():
		intersect.append("Intergenic")
		pos = v_row["POS"]
		for j, o_row in orf.iterrows():
			if o_row["Start"] <= pos and pos <= o_row["End"]:
				if intersect[i] == "Intergenic":
					intersect[i] = o_row["Gene"]
				else:
					intersect[i] += f",{o_row['Gene']}"
	out = vcf.rowdata[["CHROM", "POS", "REF", "ALT"]].copy()
	out.rename(columns={"CHROM": "Chromosome", "POS": "Position", "REF": "Reference", "ALT": "Alternative"}, inplace=True)
	out["ORF"] = intersect
	return out


def annotate(fasta, out_dir, gbk_fn, ref_fn, snpeff, verbose, internal):

	print("##################")
	print("# Annotate FASTA #")
	print("##################")
	print(f"FASTA: {fasta}")
	print(f"Output: {out_dir}")
	print(f"Genbank: {gbk_fn}")
	print(f"Reference: {ref_fn}")
	print(f"Verbose: {verbose}")

	qries = read_fasta(fasta)
	ref_gbk = read_ref_genbank(gbk_fn)
	for qry in qries:
		if isinstance(qry, Error) and qry.id == 0:
			os.makedirs(f"{out_dir}/{qry.qryid}/", exist_ok=True)
			with open(f"{out_dir}/{qry.qryid}/{qry.qryid}.log", "w") as f:
				f.write(qry.msg + "\n")
			continue
		anno = Annotation(qry, ref_gbk, out_dir, verbose)
		anno.run()
		# anno.anno_df.to_csv(f"{out_dir}/{qry.id}/{qry.id}.tsv", sep="\t", index=False)
		anno.orf.to_csv(f"{out_dir}/{qry.id}/{qry.id}_orf.tsv", sep="\t", index=False)
		anno.write_sequences(f"{out_dir}/{qry.id}/{qry.id}")

		print("Making VCF.")
		fasta2vcf(fasta_fn=None, ref_fn=ref_fn, \
			align_fn=anno.align_fn, out_dir=f"{out_dir}/{qry.id}", \
			compress_vcf=False, clean_up=True, verbose=False)

		if internal:
			intersect = vcf_intersect_orf(f"{out_dir}/{qry.id}/{qry.id}.vcf", anno.orf)
			intersect.to_csv(f"{out_dir}/{qry.id}/{qry.id}.display.tsv", index=False, sep="\t")

		if snpeff:
			run_snpEff(f"{out_dir}/{qry.id}/{qry.id}.vcf", f"{out_dir}/{qry.id}/{qry.id}.snpEff.vcf")
			
			ORF1a_start = int(anno.orf.loc[anno.orf["Gene"] == "ORF1a","Start"].values[0])
			ORF1a_end = int(anno.orf.loc[anno.orf["Gene"] == "ORF1a","End"].values[0])
			snpEff = parse_snpEff(f"{out_dir}/{qry.id}/{qry.id}.snpEff.vcf", ORF1a=[ORF1a_start,ORF1a_end])
			snpEff.to_csv(f"{out_dir}/{qry.id}/{qry.id}.snpEff.tsv", index=False, sep="\t")


@click.command()
@click.option("-f", "--fasta", type=str, required=True, \
	help="Fasta file containing one or more virus strains.")
@click.option("-o", "--out_dir", type=str, required=False, \
	help="Output directory", default="results", show_default=True)
@click.option("--gbk_fn", "-g", type=str, required=False, \
	help="Genbank file.", default="data/NC_045512.2.gbk", \
	show_default=True)
@click.option("--ref_fn", "-r", type=str, required=False,\
	help="Reference FASTA file.", default="data/NC_045512.2.fasta", \
	show_default=True)
@click.option("--snpeff", type=bool, default=True, \
	help="Whether to run snpEff", show_default=True)
@click.option("--verbose", "-v", is_flag=True, default=False, \
	help="Verbosity")
@click.option("--internal", is_flag=True, default=False, \
	help="Internal Use.")
def main(fasta, out_dir, gbk_fn, ref_fn, snpeff, verbose, internal):
	annotate(fasta, out_dir, gbk_fn, ref_fn, snpeff, verbose, internal)


if __name__ == "__main__":
	main()
