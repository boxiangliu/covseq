import os
import argparse
from collections import defaultdict, OrderedDict
import numpy as np
import pandas as pd
from Bio import SeqIO, SeqFeature, Seq
import ipdb
parser = argparse.ArgumentParser(description='Parse virus genome annotation.')
parser.add_argument("--in_prefix", type=str,
					help="Prefix for results from VAPID.")
parser.add_argument("--out_dir", type=str,
					help="Output directory.")
args = parser.parse_args()
in_prefix = args.in_prefix
out_dir = args.out_dir

in_prefix = "/mnt/scratch/boxiang/projects/viraviz/processed_data/parse_vapid/vapid/EPI_ISL_402131/EPI_ISL_402131"
out_dir = "../processed_data/parse_vapid/parsed/"

class Annotation():
	def __init__(self, prefix):
		self.align_fn = f"{prefix}.ali"
		self.ref_fn = f"{prefix}_ref.fasta"
		self.gbf_fn = f"{prefix}.gbf"
		self.get_ref_id()
		self.get_qry_id()
		self.parse_align()
		self.get_mutation()
		self.get_coverage_and_percent_identity()
		self.parse_genbank()


	def get_ref_id(self):
		with open(self.ref_fn, "r") as f:
			for line in f:
				self.ref_id = line.strip().replace(">", "")
				break

	def get_qry_id(self):
		with open(self.align_fn, "r") as f:
			for line in f:
				self.qry_id = line.strip().replace(">", "")
				break

	def parse_align(self):
		align = defaultdict(list)
		coord = defaultdict(list)
		with open(self.align_fn, "r") as f:
			for line in f:
				if line.startswith(">"):
					strain = line.strip().replace(">", "")

					if strain == self.ref_id:
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

		self.nt_df = pd.DataFrame({
			"ref_coord": ref_coord,
			"ref_nt": align["ref"],
			"qry_coord": qry_coord,
			"qry_nt": align["qry"]
			})


	def get_mutation(self):
		# def print_mut(x, y):
		# 	if x == y:
		# 		return "same"
		# 	else:
		# 		return f"{x},{y}"

		ref_nt = self.nt_df["ref_nt"]
		qry_nt = self.nt_df["qry_nt"]
		# self.nt_df["nt_mut"] = self.nt_df.apply(\
		# 	lambda x: print_mut(x['ref_nt'], x['qry_nt']), axis=1)

		self.nt_df["nt_mut"] = [f"{i},{j}" if i!=j else "same" \
			for i,j in zip(ref_nt, qry_nt)]


	def get_coverage_and_percent_identity(self):

		nt = self.nt_df

		qry_nt = nt[nt["qry_nt"] != "-"]
		qry_seq_length = qry_nt.shape[0]
		qry_cover_length = qry_nt[qry_nt["ref_nt"] != "-"].shape[0]
		coverage = qry_cover_length / qry_seq_length
		self.coverage = coverage

		identical = (qry_nt["nt_mut"] == "same").sum()
		pct_iden = identical / qry_seq_length
		self.percent_identity = pct_iden

		ref_nt = nt[nt["ref_nt"] != "-"]
		ref_seq_length = ref_nt.shape[0]
		ref_cover_length = ref_nt[ref_nt["qry_nt"] != "-"].shape[0]
		inverse_coverage = ref_cover_length / ref_seq_length
		self.inverse_coverage = inverse_coverage


	def parse_genbank(self):

		if self.inverse_coverage < 0.7:
			print(f"Inverse coverage (={self.inverse_coverage}) is too low.")
			print(f"Unable to parse genbank information.")
			return 

		nt_df = self.nt_df

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


		def qry_location2ref_location(qry_location, nt_df):

			def get_ref_location(qry_start, qry_end):

				ref_start = int(nt_df.loc[nt_df["qry_coord"] == qry_start, "ref_coord"].values[0])
				ref_end = int(nt_df.loc[nt_df["qry_coord"] == qry_end, "ref_coord"].values[0])
				ref_location = SeqFeature.FeatureLocation(ref_start, ref_end)
			
				return ref_location

			if isinstance(qry_location, SeqFeature.FeatureLocation):
				
				ref_location = get_ref_location(qry_location.start, qry_location.end)

			elif isinstance(qry_location, SeqFeature.CompoundLocation):

				loc_0 = qry_location.__dict__["parts"][0]
				qry_start = loc_0.start
				qry_end = loc_0.end
				
				ref_location = get_ref_location(qry_start, qry_end)

				for loc_i in qry_location.__dict__["parts"][1:]:
					qry_start = loc_i.start
					qry_end = loc_i.end
					ref_location += get_ref_location(qry_start, qry_end)

			else:
				raise Exception(type(qry_location) + " not defined!")
			
			return ref_location


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

		gb = next(SeqIO.parse(self.gbf_fn, "genbank"))
		protein_container = OrderedDict()

		for feature in gb.features:
			if feature.type == "source":
				continue

			elif feature.type == "CDS":

				qry_protein_anno = get_protein_anno(feature.location, gb.seq)

				ref_location = qry_location2ref_location(feature.location, self.nt_df)
				ref_seq = Seq.Seq("".join(nt_df["ref_nt"].tolist()).replace("-",""))
				ref_protein_anno = get_protein_anno(ref_location, ref_seq)

			else:
				raise Exception(feature.type + " not defined!")

			protein_name = feature.qualifiers["product"][0].replace(" ","_")

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

		merged = self.nt_df

		for k, v in protein_container.items():
			merged = pd.merge(merged, v["ref"], on="ref_coord", how="outer")
			merged = pd.merge(merged, v["qry"], on="qry_coord", how="outer")
			ref_aa = merged[f"{k}:ref_aa"].fillna("").tolist()
			qry_aa = merged[f"{k}:qry_aa"].fillna("").tolist()

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

			merged[f"{k}:aa_mut"] = [compare_aa(i,j) \
				for i, j in zip(ref_aa, qry_aa)]

		self.anno_df = merged.sort_values("ref_coord")


if __name__ == "__main__":
	print(in_prefix)
	anno = Annotation(in_prefix)
	anno.anno_df.to_csv(f"{out_dir}/{anno.qry_id}.tsv", sep="\t", index=False)