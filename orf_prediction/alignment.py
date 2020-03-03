import sys
import pandas as pd
import os
import glob

in_dir = sys.argv[1]
out_dir = sys.argv[2]
os.makedirs(out_dir, exist_ok=True)

# fn = "../data/blast_aligned_wuhan_vs_90seqs/res_align_wuhan_vs_EPI_ISL_408670"

class Alignment():
	def __init__(self, fn):
		self.fn = fn
		self.strain = os.path.basename(fn).split("_vs_")[1]
		self.length = None
		self.alignment = []
		self.map = pd.DataFrame()
		
		self.parse()
		self.create_map()

	def parse(self):
		with open(self.fn, "r") as f:
			while True:
				line = f.readline()
				if line == "":
					break 

				elif line.startswith("Length="):
					self.length = int(line.replace("Length=", ""))

				elif line.startswith("Query  "):
					split_line = line.strip().split()
					query = {}
					query["start"] = int(split_line[1])
					query["end"] = int(split_line[3])
					query["seq"] = split_line[2]
					
					f.readline() # skip |'s
					line = f.readline()
					split_line = line.strip().split()
					sbjct = {}
					sbjct["start"] = int(split_line[1])
					sbjct["end"] = int(split_line[3])
					sbjct["seq"] = split_line[2]

					self.alignment.append({"query": query, "sbjct": sbjct})

				else:
					pass

	def create_map(self):
		query_pos = []
		sbjct_pos = []
		for alignment in self.alignment:
			query = alignment["query"]
			sbjct = alignment["sbjct"]
			query_pointer = query["start"]
			sbjct_pointer = sbjct["start"]

			query_pos.append(query_pointer)
			sbjct_pos.append(sbjct_pointer)

			for bi, bj in zip(query["seq"][1:], sbjct["seq"][1:]):
				if bi != "-":
					query_pointer += 1
					query_pos.append(query_pointer)
				else:
					query_pos.append(float("nan"))

				if bj != "-":
					sbjct_pointer += 1
					sbjct_pos.append(sbjct_pointer)
				else:
					sbjct_pos.append(float("nan"))
		self.map = pd.DataFrame({"query": query_pos, "sbjct": sbjct_pos})
		self.map["strain"] = self.strain

if __name__ == "__main__":
	container = []
	for fn in glob.glob(f"{in_dir}/res_align*"):
		alignment = Alignment(fn)
		container.append(alignment.map)
	map_ = pd.concat(container)
	map_.to_csv(f"{out_dir}/map.csv", index=False, sep="\t")