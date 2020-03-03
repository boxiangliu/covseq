import pandas as pd
import os
import ipdb

alignment_fn = "../processed_data/alignment/map.csv"
predict_fn = "../processed_data/glimmer/gisaid/iterated.predict"
out_dir = "../processed_data/coord2ref/"
os.makedirs(out_dir, exist_ok=True)

alignment = pd.read_table(alignment_fn)
alignment.set_index(["strain", "sbjct"])

with open(predict_fn, "r") as fin, \
	open(f"{out_dir}/iterated.predict.ref_coord", "w") as fout:
	for line in fin:
		if line.startswith(">"):
			strain = line.strip().replace(">","")
			print(strain)
			fout.write(line)
		else:
			split_line = line.strip().split()
			for i in [1, 2]:
				coord = int(split_line[i])
				print(coord)
				ref_coord = alignment[(alignment["strain"]==strain) & \
					(alignment["sbjct"]==coord)]["query"].values[0]
				split_line[i] = str(int(ref_coord))
			out_str = "\t".join(split_line)
			fout.write(out_str + "\n")