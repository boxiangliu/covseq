fasta_fn=$1
new_fasta_fn=$2
out_fn=$3

out_dir=$(dirname $out_fn)
mkdir -p $out_dir
cat $fasta_fn $new_fasta_fn > $out_fn