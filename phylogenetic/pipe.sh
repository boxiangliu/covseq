# Filter out sequences <99% sequence identity from the reference seq:
python3 phylogenetic/filter_distant_seq.py


# Sample sequences according to location and time:
python3 phylogenetic/sample_seq.py


# Add metadata (location and time) to the header of the sequences:
python3 phylogenetic/add_meta_to_header.py 
python3 phylogenetic/add_meta_to_header.py

# Runtime analysis:
python3 phylogenetic/runtime_vs_num_strains.py