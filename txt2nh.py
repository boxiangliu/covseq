#!/usr/bin/env python3
import sys 
in_fn = sys.argv[1]
out_fn = sys.argv[2]

with open(in_fn, "r") as fin, open(out_fn, "w") as fout:
    for line in fin:
        split_line = line.strip().split("\t")
        if split_line[0] == "tax_id":
            continue # skip header line
        symbol = split_line[5]
        description = split_line[7]
        start = int(split_line[12])
        end = int(split_line[13])
        orientation = split_line[14]
        orientation = "+" if orientation == "plus" else "-"
        frame = str(start % 3)
        frame = f"{orientation}{frame}"
        fout.write(f"{symbol}\t{start}\t{end}\t{frame}\t{description}\n")
