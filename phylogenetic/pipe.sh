# Filter out sequences <99% sequence identity from the reference seq:
python3 phylogenetic/filter_distant_seq.py


# Sample sequences according to location and time:
python3 phylogenetic/sample_seq.py


# Add metadata (location and time) to the header of the sequences:
python3 phylogenetic/add_meta_to_header.py --meta_fn "../data/gisaid/metadata/metadata.tsv" --fasta_fn "../processed_data/phylogenetic/filter_distant_seq/keep.fasta" --out_fn "../processed_data/phylogenetic/add_meta_to_header/keep.fasta"
python3 phylogenetic/add_meta_to_header.py --meta_fn "../data/gisaid/metadata/metadata.tsv" --fasta_fn "../processed_data/phylogenetic/sample_seq/sample.fasta" --out_fn "../processed_data/phylogenetic/add_meta_to_header/sample.fasta"


# Append new sequence:
bash phylogenetic/append_new_seq.sh "../processed_data/phylogenetic/add_meta_to_header/keep.fasta" "../data/BJ_CDC/4seq.fasta" "../processed_data/phylogenetic/append_new_seq/keep_and_new.fasta"
bash phylogenetic/append_new_seq.sh "../processed_data/phylogenetic/add_meta_to_header/sample.fasta" "../data/BJ_CDC/4seq.fasta" "../processed_data/phylogenetic/append_new_seq/sample_and_new.fasta"


# Construct tree:
python3 phylogenetic/construct_tree.py --in_fn ../processed_data/phylogenetic/append_new_seq/sample_and_new.fasta --out_prefix ../processed_data/phylogenetic/construct_tree/sample_and_new

# Runtime analysis:
python3 phylogenetic/runtime_vs_num_strains.py