import sys

for line in sys.stdin:
	if line.startswith("#"):
		sys.stdout.write(line)
		continue
	split_line = line.split("\t")
	if not split_line[3].endswith("AAAAAAAAAAAAAAAAAAAAAAAAAAAA"):
		sys.stdout.write(line)

