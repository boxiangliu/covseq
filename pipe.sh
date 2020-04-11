# Download COVID-19 sequences from GISAID:
# Accessed March 24, 2020
python3 download/download_gisaid.py -u lbxjollier -p 71RwYNz4nljy --out_dir /Users/boxiang/Documents/work/Baidu/projects/viraviz/data/gisaid/

# Concatenate all FASTA records:
python3 preprocess/concatenate_fasta.py
python3 preprocess/filter_fasta.py

# Profile annotation.py
bash annotation/profile_annotation.sh 
# Figure:../processed_data/annotation/profile_annotation/stat.pdf
# Result: The merge operation in transfer_feature is the slowest, consuming about 33% of overall runtime
# Conclusion: 33% improve does not justify spending time. Leaving annotation.py as it is.

# Create polygenetic tree:


/mnt/scratch/boxiang/projects/open_reading_frame/src/glimmer3.02/scripts/g3-from-training.csh \
    /mnt/scratch/boxiang/projects/open_reading_frame/data/COVID19_dna.txt \
    /mnt/scratch/boxiang/projects/open_reading_frame/data/COVID19.nh \
    ../processed_data/glimmer/COVID19.from-training

/mnt/scratch/boxiang/projects/open_reading_frame/src/glimmer3.02/scripts/g3-iterated.csh \
    /mnt/scratch/boxiang/projects/open_reading_frame/data/COVID19_dna.txt \
    ../processed_data/glimmer/COVID19.iterated


/mnt/scratch/boxiang/projects/open_reading_frame/src/glimmer3.02/scripts/g3-iterated.csh \
    ../processed_data/prepare_g3_input/g3_input.fasta \
    ../processed_data/glimmer/gisaid/iterated


python3 alignment.py ../data/blast_aligned_wuhan_vs_90seqs ../processed_data/alignment/

python3 coord2ref.py


# Converting FASTA to VCF: 
python3 fasta2vcf/fasta2vcf.py


# Plot the distribution of mutations:
python3 fasta2vcf/plot_mut_per_sample.py
# Figure: ../processed_data/fasta2vcf/mut_per_sample/mut_per_sample.pdf
# Figure: ../processed_data/fasta2vcf/mut_per_sample/mut_per_sample_lt_200.pdf
# Conclusion: remove samples with > 50 mutations.


# Filter samples with > 50 mutations:
python3 fasta2vcf/filter_samples.py


# Merge VCF files:
bash merge_vcfs/merge_vcfs.sh

# Filter sites with >2 alleles
# Also filter sites within the poly-A tail
bash merge_vcfs/filter_sites.sh

# Plot the mutation distribution:
python3 mutation_distribution/plot_sample_per_mutation.py

# Plot mutation over time:
python3 mutation_distribution/plot_mutation_over_time.py

# Make phenotype table:
python3 phenotype/phenotype.py