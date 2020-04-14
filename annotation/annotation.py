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


def read_fasta(fasta_fn):
	fasta = list(SeqIO.parse(fasta_fn, "fasta"))
	for x in fasta:
		x.id = x.id.replace(" ", "_").replace("/", "_")
	return fasta


def read_ref_genbank():
	return next(SeqIO.parse("data/ref.gbk", "genbank"))


class Annotation():
	def __init__(self, qry, ref_gbk, out_dir):
		'''
		qry: a SeqIO FASTA object.
		ref_gbk: a SeqIO genbank object.
		out_dir: output directory.
		'''
		self.qry = qry
		self.ref_gbk = ref_gbk
		self.out_dir = out_dir

	def run(self):
		print("Start annotating the query sequence...")
		print(f"Query ID: {self.qry.id}")
		print(f"Ref. ID: {self.ref_gbk.id}")
		if len(self.qry.seq) < 29000:
			print("Query sequence should have at least 29000 nucleotides! Aborting...")
			return 

		self.align(self.ref_gbk, self.qry, self.out_dir)
		self.parse_align(self.align_fn)
		self.get_mutation(self.nt_df)
		self.transfer_feature(self.ref_gbk, self.nt_df)


	def align(self, ref, qry, out_dir):
		print("Aligning reference and query sequences...")
		work_dir = f"{out_dir}/{qry.id}"
		os.makedirs(work_dir, exist_ok=True)
		
		op_sys = subprocess.run("uname", capture_output=True).stdout
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
		output = subprocess.run(cmd, shell=True, capture_output=True)

		self.align_fn = mafft_out_fn


	def parse_align(self, align_fn):
		print("Parsing alignment...")
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

					counter = 0
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
		print("Finding nucleotide mutations between ref. and query sequences...")
		ref_nt = nt_df["ref_nt"]
		qry_nt = nt_df["qry_nt"]

		nt_df["nt_mut"] = [f"{i},{j}" if i!=j else "same" \
			for i,j in zip(ref_nt, qry_nt)]

		self.nt_df = nt_df


	def transfer_feature(self, ref_gbk, nt_df):
		print("Starting transfering CDS fetaures from ref to query...")
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


		def get_protein_anno(location, seq):
			protein_anno = init_protein_annotation(location)

			if isinstance(location, SeqFeature.FeatureLocation):
				
				protein_anno = update_protein_annotation(location, seq, protein_anno)

			elif isinstance(location, SeqFeature.CompoundLocation):
				
				for location in location.__dict__["parts"]:
					protein_anno = update_protein_annotation(location, seq, protein_anno)

			else: 

				raise Exception(type(location) + " not defined!")

			return protein_anno


		def transfer_simple_location(start, end, fro, to):
			new_start = int(to[fro == start].values[0])
			new_end = int(to[fro == end].values[0])
			location = SeqFeature.FeatureLocation(new_start, new_end)
			return location


		def transfer_location(location, fro, to):

			if isinstance(location, SeqFeature.FeatureLocation):
				# FeatureLocation has only 1 start and 1 end.
				new_location = transfer_simple_location(location.start, \
					location.end, fro, to)

			elif isinstance(location, SeqFeature.CompoundLocation):
				# CompoundLocation has multiple locations.
				loc_0 = location.__dict__["parts"][0]
				new_location = transfer_simple_location(loc_0.start, \
					loc_0.end, fro, to)

				for loc_i in location.__dict__["parts"][1:]:
					new_location += transfer_simple_location(loc_i.start, \
						loc_i.end, fro, to)

			else:
				raise Exception(type(location) + " not defined!")
			
			return new_location


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


		protein_container = OrderedDict()

		for ref_feature in ref_gbk.features:

			if ref_feature.type == "CDS":

				protein_name = ref_feature.qualifiers["product"][0].replace(" ","_")
				print(f" - ORF: {protein_name}")
				
				print(f" - Annotating reference protein...")
				ref_protein_anno = get_protein_anno(ref_feature.location, ref_gbk.seq)

				print(f" - Transfering features from reference to query...")
				qry_location = transfer_location(ref_feature.location, \
					nt_df["ref_coord"], nt_df["qry_coord"])
				qry_seq = Seq.Seq("".join(nt_df["qry_nt"].tolist()).replace("-",""))
				qry_protein_anno = get_protein_anno(qry_location, qry_seq)


				print(" - Converting protein annotations into DataFrames...")
				colnames = ["qry_coord", f"{protein_name}:qry_aa",f"{protein_name}:frame",
					f"{protein_name}:ribo_slip",f"{protein_name}:rna_edit"]
				qry_protein_df = dict2df(qry_protein_anno, colnames)
				qry_protein_df["qry_coord"] = qry_protein_df["qry_coord"].astype("Int64")


				colnames = ["ref_coord", f"{protein_name}:ref_aa", "placeholder1", 
					"placeholder2", "placeholder3"]
				ref_protein_df = dict2df(ref_protein_anno, colnames)
				ref_protein_df["ref_coord"] = ref_protein_df["ref_coord"].astype("Int64")
				ref_protein_df.drop(["placeholder1", "placeholder2", "placeholder3"], axis=1, inplace=True)

				protein_container[protein_name]={"ref": ref_protein_df, "qry": qry_protein_df}


		print("Finding amino acid mutations...")
		merged = nt_df

		for k, v in protein_container.items():
			merged = pd.merge(merged, v["ref"], on="ref_coord", how="outer")
			merged = pd.merge(merged, v["qry"], on="qry_coord", how="outer")
			ref_aa = merged[f"{k}:ref_aa"].fillna("").tolist()
			qry_aa = merged[f"{k}:qry_aa"].fillna("").tolist()

			merged[f"{k}:aa_mut"] = [compare_aa(i,j) \
				for i, j in zip(ref_aa, qry_aa)]

		self.anno_df = merged.sort_values("ref_coord")


@click.command()
@click.option("-f", "--fasta", type=str, required=True, \
	help="Fasta file containing one or more virus strains.")
@click.option("-o", "--out_dir", type=str, required=False, \
	help="Output directory", default="results", show_default=True)
def annotate(fasta, out_dir):
	qries = read_fasta(fasta)
	ref_gbk = read_ref_genbank()
	for qry in qries:
		anno = Annotation(qry, ref_gbk, out_dir)
		anno.run()
		anno.anno_df.to_csv(f"{out_dir}/{qry.id}/{qry.id}.tsv", sep="\t", index=False)
		print(anno.align_fn)
		fasta2vcf(fasta_fn=None, ref_fn=None, \
			align_fn=anno.align_fn, out_dir=f"{out_dir}/{qry.id}", \
			compress_vcf=False, clean_up=True, verbose=False)

if __name__ == "__main__":
	annotate()
