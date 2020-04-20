# CoV-Seq: COVID-19 Genomic Sequence Database and Visualization

<div style="text-align:center"><img src="website/figure/1x/logo.png"/></div>

## Clone the repository
git clone https://github.com/boxiangliu/covseq.git

CoV-Seq is written and tested in `python 3.7`, and should in theory work with any `python 3` versions. CoV-Seq is not compatible with `python 2`.

## Quick Start

We have made CoV-Seq very easy to use. To annotate a SARS-CoV-2 genome in FASTA format, go to the `covseq` directory and run: 

```
python annotation/annotation.py --fasta data/GZ8H0001.fasta --out_dir results
```

Replace `data/GZ8H0001.fasta` with your fasta file and `results` with your desired output directory. 

If the command was successful, there should 5 files in the `results/` directory.

```
Guangzhou_GZ8H0001_2020.tsv
Guangzhou_GZ8H0001_2020_orf.tsv
Guangzhou_GZ8H0001_2020.vcf
Guangzhou_GZ8H0001_2020.snpEff.vcf
Guangzhou_GZ8H0001_2020.snpEff.tsv
```

The filename `Guangzhou_GZ8H0001_2020` comes from the FASTA header in `data/GZ8H0001.fasta`. Notice that characters "/" have been automatically replaced with "\_" to save us from using escape characters. Here is a description of each file:

1. \*\_orf.tsv: ORF predictions.
2. \*.vcf: variant calls
3. \*.snpEff.vcf: variant calls annotated with snpEff
4. \*.snpEff.tsv: parsed VCF annotations

**That's it**. You can now use the VCF and annotations for downstream analysis. 

For all options, run 

```
annotation/annotation.py --help
```

## Replicating CoV-Seq Results

To replicate all results in on `covseq.baidu.com/browse`, please follow instructions below. 


### Preliminaries 

You will need to install `bcftools` and `htslib` by following instructions [here](http://www.htslib.org/download/). 


### Download genomic data and metadata

The first step is to download data from repositories. All steps assume `covseq` repo as the working directory. 

To download sequence and metadata from GISAID: 

1. Go to https://www.gisaid.org/
2. Click on login on the upper right corner
3. To download FASTA sequences: Click on EpiCoV -> Browse -> Download (lower right corner) -> Move the FASTA sequences to `../data/gisaid/fasta/`
4. Download the Acknowledgement Table (link beneath the virus display table) -> Move the table to `../data/gisaid/metadata/acknowledgement/`


To download data from NCBI: 

1. Go to https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049
2. Click "Select Columns" (upper-right corner) -> Add all columns -> Apply
3. To download FASTA sequences: Download -> Select Nucleotide (under Sequence data) -> Select Download All Records -> Select Use default -> Download -> Move FASTA files into `../data/ncbi/fasta/`
4. To download metadata: Download -> Select CSV format (under Current table view result) -> Select Download All Records -> Click Select All -> Download -> Move metadata into `../data/ncbi/metadata/`


To download data from EMBL:

1. Go to https://www.ebi.ac.uk/ena/pathogens/covid-19
2. Click Assembled Sequences (under Summary)
3. Click FASTA to download sequence data -> Move FASTA into `../data/embl/fasta/`
4. Click EMBL to download metadata -> Move metadata into `../data/embl/metadata/`


To download data from CNGB: 

1. Go to https://db.cngb.org/virus/ncov
2. Click on Data source platform once to rank CNGB first
3. Select CNGB samples 
4. Click on Download Sequence data to download FASTA sequences -> Move FASTA sequences to `../data/cngb/fasta/`
5. Click on Download Excel to download metadata -> Move metadata to `../data/cngb/metadata/`


### Preprocess

We will preprocess the data by first concatenating all FASTA files and standardize FASTA headers 

1. Concatenate all FASTA files standardize FASTA headers

```
python3 preprocess/concatenate_fasta.py -i ../data/ -o ../data/aggregated/fasta/raw.fasta --cgnb_metadata <../data/cngb/metadata/CNGBdb_VirusDIP.csv>
```

Replace `../data/cngb/metadata/CNGBdb_VirusDIP.csv` with your own CNGB metadata file.

Next we will remove incomplete genomes (number of nucleotide < 25000). 

2. Filter out incomplete genomes. 

```
python3 preprocess/filter_fasta.py -i ../data/aggregated/fasta/raw.fasta --out_dir ../processed_data/preprocess/filter_fasta/ --final_fn ../data/aggregated/fasta/preprocessed.fasta
```

Note that this command will create a `../processed_data/preprocess/filter_fasta/` to store intermedite files. 


### Call variants

Now we have preprocessed FASTA files, let's call variants from nucleotide sequences, using the RefSeq sequence `NC_045512.2` as the reference.

1. Call variants 
```
python3 vcf/fasta2vcf.py -f ../data/aggregated/fasta/preprocessed.fasta -r data/NC_045512.2.fasta -o ../data/aggregated/vcf/individual/
```

This command will create a directory called `../data/aggregated/vcf/individual/`, where VCF files from individual FASTA files will reside.


Some FASTA files have large numbers of sequencing errors and will produce abnormally long lists of variants. Let's filter them out.

2. Filter out variants with > 150 mutations. 
```
python3 vcf/filter_samples.py -i ../data/aggregated/vcf/individual/ -o ../processed_data/vcf/filter_samples/ -c 150
```

Here we have chosen 150 as a default because a clear gap exists between samples with <150 mutations and samples with >1000 mutations with nothing in between. However, feel free to adjust this parameter to fit your needs. 


Next we need to merge VCF files and normalize each record. 

3. Merge VCF files

```
bash vcf/merge_vcfs.sh ../processed_data/vcf/filter_samples/ ../data/aggregated/vcf/merged/ data/NC_045512.2.fasta
```

Finally we will remove multi-allelic variants and variants within the poly-A tail. 

4. Filter variants

```
bash vcf/filter_sites.sh ../data/aggregated/vcf/merged/merged.vcf.gz ../data/aggregated/vcf/merged/filtered
```




### Annotate VCF files 



### Merging metadata 



## Frequently Asked Questions
Please see [FAQ](example.com) here. 

## Bug report 

We welcome bug report [here](https://github.com/boxiangliu/covseq/issues). Please help us by providing as much information as possible.

If you are running the source code, please provide the following:

- A minimum example to reproduce the bug;
- Expected result;
- Full error message;
- Your operating system.


If you are using the CoV-Seq website, please provide the following:

- Your input/uploads to the website;
- Step-by-step instructions to reproduce the error;
- A screenshot of the error message;
- Your browser information.