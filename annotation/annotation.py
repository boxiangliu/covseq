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
import platform 
import shutil
import json

class Status():
	'''
	Status Code:
	0: success
	1: sequence contains more than 5% non-A[T/U]GC letters
	2: sequence contains < 25,000 nucleotides or > 35000 nucleotides
	3: mafft not found
	4: unknown error
	5: empty FASTA file 
	6: empty sequence (e.g. >no_sequence\n)
	'''
	def __init__(self, code, msg, id):
		self.code = code
		self.msg = msg
		self.id = id


def read_fasta(fasta_fn):
	fasta = list(SeqIO.parse(fasta_fn, "fasta"))
	for i, x in enumerate(fasta):
		x.id = secure_filename(x.id)
		print("Replacing U's with T's.")
		x.seq = Seq.Seq(str(x.seq).upper().replace("U","T"))


		if len(x.seq) == 0:
			err_msg = f"ERROR: Sequence {x.id} is empty."
			print(err_msg)
			fasta[i] = Status(6, err_msg, x.id)
			continue 

		non_ATGC = 0
		for nuc in x.seq:
			if nuc not in "ATGC":
				non_ATGC += 1
		
		if non_ATGC/len(x.seq) > 0.05:
			err_msg = f"ERROR: Sequence {x.id} contains more than 5% non-A[T/U]GC letters."
			print(err_msg)
			fasta[i] = Status(1, err_msg, x.id)

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
		if len(self.qry.seq) < 25000 or len(self.qry.seq) > 35000:
			err_msg = f"ERROR: Sequence {self.qry.id} should be between 25,000 and 35,000 nucleotides!"
			print(err_msg)
			return Status(2, err_msg, self.qry.id)

		align_status = self.align(self.ref_gbk, self.qry, self.out_dir)
		if align_status.code == 3:
			return align_status
		self.parse_align(self.align_fn)
		# self.get_mutation(self.nt_df)
		# self.transfer_feature(self.ref_gbk, self.nt_df, self.verbose)
		self.get_orf(self.ref_gbk, self.nt_df)
		return Status(0, "success", self.qry.id)

	def align(self, ref, qry, out_dir):
		print("Aligning reference and query sequences.")
		work_dir = f"{out_dir}/{qry.id}"
		os.makedirs(work_dir, exist_ok=True)
		
		if platform.system() == "Linux":
			mafft = "ext/mafft-linux64/mafft.bat"
		elif platform.system() == "Darwin":
			mafft = "ext/mafft-mac/mafft.bat"
		elif platform.system() == "Windows":

			cmd = 'mafft'
			whereis_mafft = subprocess.check_output('whereis '+cmd, stderr=subprocess.STDOUT, shell=True).decode("utf-8").rstrip()
			if whereis_mafft == cmd+':':
				err_msg = "ERROR: mafft not found! Install it here: https://mafft.cbrc.jp/alignment/software/windows.html"
				return Status(3, err_msg, self.qry.id)
			else:
				mafft = 'bash '+cmd

		mafft_in_fn = f"{work_dir}/{qry.id}.fasta"
		with open(mafft_in_fn, "w") as f:
			f.write(f">{ref.id}\n")
			f.write(f"{str(ref.seq)}\n")
			f.write(f">{qry.id}\n")
			f.write(f"{str(qry.seq)}\n")

		mafft_out_fn = f"{work_dir}/{qry.id}.ali"
		cmd = f"{mafft} {mafft_in_fn} > {mafft_out_fn}"
		self.align_fn = mafft_out_fn

		if platform.system() == "Linux" or platform.system() == "Darwin":
			output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True).decode("utf-8")

			if "Error" in output or "Warning" in output:
				err_msg = f"Sequence {self.qry.id}: {output}"
				return Status(3, err_msg, self.qry.id)
			else:
				return Status(0, "success", self.qry.id)

		elif platform.system() == "Windows":
			output = os.system(cmd)
			if output == 0:
				return Status(0, "success", self.qry.id)
			elif output == 1:
				err_msg = f"Unknown error in Windows OS for {self.qry.id}: {output}" # TODO: may change the message here
				return Status(4, err_msg, self.qry.id)


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
		container = defaultdict(dict)

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
			container["nt"][gene] = str(nt_seq)

			if "," in start:
				aa_starts = [int(x) for x in start.split(",")]
				aa_ends = [int(x) for x in end.split(",")]
				aa_seq = ""
				for aa_start, aa_end in zip(aa_starts, aa_ends):
					nt_seq = qry.seq[(aa_start-1):aa_end]
					aa_seq += nt_seq.translate()
			else:
				aa_seq = nt_seq.translate()
			container["aa"][gene] = str(aa_seq)

		with open(f"{out_prefix}_seq.json", "w") as f:
			json.dump(container, f)

def run_snpEff(vcf_fn, out_fn):
	print("Running snpEff.")
	cmd = f"java -jar ext/snpEff/snpEff.jar NC_045512.2 {vcf_fn} > {out_fn}"

	if platform.system() == "Linux" or platform.system() == "Darwin":
		output = subprocess.check_output(cmd, stderr=subprocess.DEVNULL, shell=True)
	elif platform.system() == "Windows":
		output = os.system(cmd)
		if output == 0:
			print("snpEff succeeded.")
			os.remove("snpEff_genes.txt")
			os.remove("snpEff_summary.html")
		else:
			print("snpEff failed. Please file an issue with your VCF file at https://github.com/boxiangliu/covseq/issues")




def vcf_intersect_orf(vcf_fn, orf):
	intersect = []
	vcf = VCF(vcf_fn)
	for i, v_row in vcf.rowdata.iterrows():
		intersect.append("Intergenic")
		pos = v_row["POS"]
		for j, o_row in orf.iterrows():
			if "," in o_row["Start"]:
				start = int(o_row["Start"].split(",")[0])
				end = int(o_row["End"].split(",")[-1])
			else:
				start = int(o_row["Start"])
				end = int(o_row["End"])

			if start <= pos and pos <= end:
				if intersect[i] == "Intergenic":
					intersect[i] = o_row["Gene"]
				else:
					intersect[i] += f",{o_row['Gene']}"
	out = vcf.rowdata[["CHROM", "POS", "REF", "ALT"]].copy()
	out.rename(columns={"CHROM": "Chromosome", "POS": "Position", "REF": "Reference", "ALT": "Alternative"}, inplace=True)
	out["ORF"] = intersect
	return out


def write_error(out_dir, status):
	if status.code == 5: # empty fasta
		with open(f"{out_dir}/{status.id}.err", "w") as f:
			f.write(status.msg + "\n")
	else:
		with open(f"{out_dir}/{status.id}/{status.id}.err", "w") as f:
			f.write(status.msg + "\n")


def annotate(fasta, out_dir, gbk_fn, ref_fn, snpeff, verbose, internal, debug):

	print("##################")
	print("# Annotate FASTA #")
	print("##################")
	print(f"FASTA: {fasta}")
	print(f"Output: {out_dir}")
	print(f"Genbank: {gbk_fn}")
	print(f"Reference: {ref_fn}")
	print(f"Verbose: {verbose}")

	if debug:
		clean_up = False
	else:
		clean_up = True

	qries = read_fasta(fasta)

	if len(qries) == 0:
		err_msg = "ERROR: FASTA file is empty!"
		empty_status = Status(5, err_msg, os.path.basename(fasta))
		write_error(out_dir, empty_status)
		print(err_msg)
		return 

	ref_gbk = read_ref_genbank(gbk_fn)
	for qry in qries:
		os.makedirs(f"{out_dir}/{qry.id}/", exist_ok=True)

		if isinstance(qry, Status) and (qry.code == 1 or qry.code == 6): # too many non-ATCG nucleotides or sequence is empty
			write_error(out_dir, qry)
			continue

		anno = Annotation(qry, ref_gbk, out_dir, verbose)
		anno_status = anno.run()

		if isinstance(anno_status, Status) and anno_status.code != 0:
			write_error(out_dir, anno_status)
			continue

		try:
			# anno.anno_df.to_csv(f"{out_dir}/{qry.id}/{qry.id}.tsv", sep="\t", index=False)
			anno.orf.to_csv(f"{out_dir}/{qry.id}/{qry.id}_orf.tsv", sep="\t", index=False)

			print("Making VCF.")
			fasta2vcf(fasta_fn=None, ref_fn=ref_fn, \
				align_fn=anno.align_fn, out_dir=f"{out_dir}/{qry.id}", \
				compress_vcf=False, clean_up=clean_up, verbose=verbose)

			if internal:
				anno.write_sequences(f"{out_dir}/{qry.id}/{qry.id}")
				intersect = vcf_intersect_orf(f"{out_dir}/{qry.id}/{qry.id}.vcf", anno.orf)
				intersect.to_csv(f"{out_dir}/{qry.id}/{qry.id}.display.tsv", index=False, sep="\t")

			if snpeff:
				run_snpEff(f"{out_dir}/{qry.id}/{qry.id}.vcf", f"{out_dir}/{qry.id}/{qry.id}.snpEff.vcf")
				try:
					ORF1a_start = int(anno.orf.loc[anno.orf["Gene"] == "ORF1a","Start"].values[0])
					ORF1a_end = int(anno.orf.loc[anno.orf["Gene"] == "ORF1a","End"].values[0])
					snpEff = parse_snpEff(f"{out_dir}/{qry.id}/{qry.id}.snpEff.vcf", ORF1a=[ORF1a_start,ORF1a_end])
				except:
					snpEff = parse_snpEff(f"{out_dir}/{qry.id}/{qry.id}.snpEff.vcf", ORF1a=None)
				snpEff.to_csv(f"{out_dir}/{qry.id}/{qry.id}.snpEff.tsv", index=False, sep="\t")
		except:
			err_msg = f"Unknown Error for sequence {qry.id}! Please submit a bug report with your input file at https://github.com/boxiangliu/covseq/issues"
			other_status = Status(4, err_msg, qry.id)
			write_error(out_dir, other_status)


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
@click.option("--debug", is_flag=True, default=False, \
	help="Activate debug mode.")
def main(fasta, out_dir, gbk_fn, ref_fn, snpeff, verbose, internal, debug):
	annotate(fasta, out_dir, gbk_fn, ref_fn, snpeff, verbose, internal, debug)


if __name__ == "__main__":
	main()
