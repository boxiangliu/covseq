fasta_fn="../processed_data/phylogenetic/add_meta_to_header/keep.fasta"
new_fasta_fn="../data/BJ_CDC/4seq.fasta"
out_dir="../processed_data/phylogenetic/append_new_seq/"
mkdir -p $out_dir
cat $fasta_fn $new_fasta_fn > $out_dir/keep_and_new.fasta