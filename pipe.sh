# Download COVID-19 sequences from GISAID:
# Accessed March 24, 2020
python3 download/download_gisaid.py -u lbxjollier -p 71RwYNz4nljy --out_dir /Users/boxiang/Documents/work/Baidu/projects/viraviz/data/gisaid/

# Concatenate all FASTA records.
# Standardize all headers.
python3 preprocess/concatenate_fasta.py -i ../data/ -o ../data/aggregated/fasta/raw.fasta --cgnb_metadata "../data/cngb/metadata/CNGBdb_VirusDIP_excel20200411_all(24)_57b4ce53c4d6c49c978596677a112211.csv"

# Filter FASTA records.
# Remove duplicates.
# Filter for complete genomes.
python3 preprocess/filter_fasta.py -i ../data/aggregated/fasta/raw.fasta --out_dir ../processed_data/preprocess/filter_fasta/ --final_fn ../data/aggregated/fasta/preprocessed.fasta


# Converting FASTA to VCF: 
python3 vcf/fasta2vcf.py -f ../data/aggregated/fasta/preprocessed.fasta -r ../data/reference/NC_045512.2.fasta -o ../data/aggregated/vcf/individual/


# Plot the distribution of mutations:
python3 vcf/filter_samples.py -i ../data/aggregated/vcf/individual/ -o ../processed_data/vcf/filter_samples/ -c 150
# Conclusion: remove samples with > 150 mutations.

# Merge VCF files:
bash vcf/merge_vcfs.sh ../processed_data/vcf/filter_samples/ ../data/aggregated/vcf/merged/ ../data/reference/NC_045512.2.fasta

# Filter sites with >2 alleles
# Also filter sites within the poly-A tail
bash vcf/filter_sites.sh ../data/aggregated/vcf/merged/merged.vcf.gz ../data/aggregated/vcf/merged/filtered

# Plot the mutation distribution:
python3 mutation_distribution/plot_sample_per_mutation.py

# Plot mutation over time:
python3 mutation_distribution/plot_mutation_over_time.py

# Make phenotype table:
python3 phenotype/parse_gisaid_metadata.py --metadata_dir ../data/gisaid/metadata/ -o ../data/aggregated/metadata/individual/gisaid.tsv
python3 phenotype/parse_ncbi_metadata.py --gb_fn ../data/ncbi/metadata/sequence.gb -o ../data/aggregated/metadata/individual/ncbi.tsv
python3 phenotype/parse_embl_metadata.py --embl_fn ../data/embl/metadata/ena_sequence_update_20200411-0207.txt -o ../data/aggregated/metadata/individual/embl.tsv
python3 phenotype/rename_cngb_metadata.py --in_fn "../data/cngb/metadata/CNGBdb_VirusDIP_excel20200411_all(24)_57b4ce53c4d6c49c978596677a112211.csv" --out_fn ../data/aggregated/metadata/individual/cngb.tsv

# Profile annotation.py
bash annotation/profile_annotation.sh 
# Figure:../processed_data/annotation/profile_annotation/stat.pdf
# Result: The merge operation in transfer_feature is the slowest, consuming about 33% of overall runtime
# Conclusion: 33% improve does not justify spending time. Leaving annotation.py as it is.
