import pandas as pd
from Bio import SeqIO
import os
import click


@click.command()
@click.option("--meta_fn", help="Metadata file.")
@click.option("--fasta_fn", help="Fasta file.")
@click.option("--out_fn", help="Output file.")
def main(meta_fn, fasta_fn, out_fn):
	out_dir = os.path.dirname(out_fn)
	os.makedirs(out_dir, exist_ok=True)

	meta = pd.read_table(meta_fn)
	meta["strain2"] = meta.apply(lambda x: x["strain"].replace(" ", "_"), axis=1)
	meta["new_header"] = meta.apply(lambda x: f'{x["strain2"]}|{x["region"]}|{x["country"]}|{x["date"]}', axis=1)
	old2new = dict()
	for row in meta.itertuples():
		old2new[row.strain] = row.new_header

	with open(out_fn, "w") as f:
		for record in SeqIO.parse(fasta_fn, "fasta"):
			record.id = old2new[record.id]
			record.description = old2new[record.description]
			SeqIO.write(record, f, "fasta")


if __name__ == "__main__":
	main()