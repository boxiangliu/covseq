# Download COVID-19 sequences from GISAID:
# Accessed March 24, 2020
# python3 download/download_gisaid.py -u <username> -p <password> --out_dir /Users/boxiang/Documents/work/Baidu/projects/viraviz/data/gisaid/metadata/

# Concatenate all FASTA records.
# Standardize all headers.
python3 preprocess/concatenate_fasta.py -i ../data/ -o ../data/aggregated/fasta/raw.fasta --cgnb_metadata ../data/cngb/metadata/CNGBdb_VirusDIP.csv --gisaid_metadata ../data/gisaid/metadata/metadata.tsv

# Filter FASTA records.
# Remove duplicates.
# Filter for complete genomes.
python3 preprocess/filter_fasta.py -i ../data/aggregated/fasta/raw.fasta --out_dir ../processed_data/preprocess/filter_fasta/ --final_fn ../data/aggregated/fasta/preprocessed.fasta
Rscript preprocess/plot_genome_lengths.R

# Converting FASTA to VCF:
python3 vcf/fasta2vcf.py -f ../data/aggregated/fasta/preprocessed.fasta -r data/NC_045512.2.fasta -o ../data/aggregated/vcf/individual/
python3 vcf/count_mutations_per_sample.py --vcf_dir ../data/aggregated/vcf/individual/ --out_dir ../processed_data/vcf/count_mutations_per_sample/

# Remove samples with too many mutations:
python3 vcf/filter_samples.py -i ../data/aggregated/vcf/individual/ -o ../processed_data/vcf/filter_samples/ -c 150
Rscript vcf/plot_mutation_count.R 
# Conclusion: remove samples with > 150 mutations.

# Merge VCF files:
bash vcf/merge_vcfs.sh ../processed_data/vcf/filter_samples/ ../data/aggregated/vcf/merged/ data/NC_045512.2.fasta
python vcf/count_mutation_per_site.py --vcf_fn ../data/aggregated/vcf/merged/merged.vcf.gz --out_dir ../processed_data/vcf/count_mutation_per_site/
Rscript vcf/plot_mutation_per_site.R 

# Filter sites with >2 alleles
# Also filter sites within the poly-A tail
bash vcf/filter_sites.sh ../data/aggregated/vcf/merged/merged.vcf.gz ../data/aggregated/vcf/merged/filtered

# Annotation mutation:
bash snpEff/snpEff.sh ../data/aggregated/vcf/merged/filtered.vcf.gz ../data/aggregated/vcf/merged/annotated
python3 snpEff/parse_snpEff.py --vcf_fn ../data/aggregated/vcf/merged/annotated.vcf.gz --out_fn ../processed_data/snpEff/parse_snpEff/annotated.tsv

# Plot the mutation distribution:
# python3 mutation_distribution/plot_sample_per_mutation.py

# Plot mutation over time:
# python3 mutation_distribution/plot_mutation_over_time.py

# Make phenotype table:
# python3 metadata/parse_gisaid_metadata.py --metadata_dir ../data/gisaid/metadata/detail/ -o ../data/aggregated/metadata/gisaid_detail.tsv --type detail
# python3 metadata/parse_gisaid_metadata.py --metadata_fn ../data/gisaid/metadata/gisaid_hcov-19.xls -o ../data/aggregated/metadata/gisaid_acknowledgement.tsv --type acknowledgement
python3 metadata/parse_gisaid_metadata.py --metadata_fn ../data/gisaid/metadata/metadata.tsv -o ../data/aggregated/metadata/gisaid.tsv --type nextmeta
# python3 metadata/parse_ncbi_metadata.py --gb_fn ../data/ncbi/metadata/sequence.gb -o ../data/aggregated/metadata/ncbi.tsv
python3 metadata/parse_ncbi_metadata.py --csv_fn ../data/ncbi/metadata/sequences.csv -o ../data/aggregated/metadata/ncbi.tsv
python3 metadata/parse_embl_metadata.py --embl_fn ../data/embl/metadata/ena_sequence.txt -o ../data/aggregated/metadata/embl.tsv
python3 metadata/rename_cngb_metadata.py --in_fn ../data/cngb/metadata/CNGBdb_VirusDIP.csv --out_fn ../data/aggregated/metadata/cngb.tsv
python3 metadata/merge_metadata.py --in_dir ../data/aggregated/metadata/ --out_prefix ../data/aggregated/metadata/merged --genome_length_fn ../processed_data/preprocess/filter_fasta/genome_lengths.tsv --duplicate_seq_fn ../processed_data/preprocess/filter_fasta/duplicate_seq.tsv --num_variant_fn ../processed_data/vcf/count_mutations_per_sample/mutations_per_sample.tsv --vcf_fn ../data/aggregated/vcf/merged/merged.vcf.gz

# Zip up VCF files and metadata files:
bash aws_upload/upload.sh


# Generate HTML from Rmarkdown:
Rscript website/rmd2html.R /Users/boxiang/Documents/work/Baidu/projects/covseq/scripts/website/faq.Rmd
Rscript website/rmd2html.R /Users/boxiang/Documents/work/Baidu/projects/covseq/scripts/website/browse.Rmd
Rscript website/rmd2html.R /Users/boxiang/Documents/work/Baidu/projects/covseq/scripts/website/index.Rmd
# MSA / Phylogenetic tree:
# python3 phylogenetic/construct_tree.py --in_fn ../data/aggregated/fasta/preprocessed.fasta --out_dir ../data/aggregated/msa/


# Profile annotation.py
# bash annotation/profile_annotation.sh 
# Figure:../processed_data/annotation/profile_annotation/stat.pdf
# Result: The merge operation in transfer_feature is the slowest, consuming about 33% of overall runtime
# Conclusion: 33% improve does not justify spending time. Leaving annotation.py as it is.
