import click
from datetime import date
import time

@click.group()
def coviz():
	pass

@coviz.command()
@click.option("-f", "--fasta", type=str, required=True, \
	help="Fasta file containing one or more virus strains.")
@click.option("-d", "--date", \
	type=click.DateTime(formats=["%Y-%m-%d"]), \
	default=str(date.today()), show_default=True, \
	help="When was the virus sequenced?")
@click.option("-r", "--region", type=str, default="Sunnyvale, CA", \
	show_default=True, help="Where was the virus obtained?")
@click.option("-c", "--coverage", type=int, default=0, \
	show_default=True, help="What is the sequencing depth?")

def annotate(fasta, date, region, coverage):
	click.echo(f"FASTA file: {fasta}")
	click.echo(f"date: {date}")
	click.echo(f"region: {region}")
	click.echo(f"coverage: {coverage}")

# @coviz.command()
# def predict():
# 	pass

# @coviz.command()
# def phylogeny():
# 	pass


if __name__ == "__main__":
	coviz()