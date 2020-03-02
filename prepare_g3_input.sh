in_dir=../data/GISAID_Human_Coronavirus_fasta_clean/
threshold=25000
out_dir=../processed_data/prepare_g3_input/
mkdir -p $out_dir

[[ -f $out_dir/g3_input.fasta ]] && rm $out_dir/g3_input.fasta
for fn in $in_dir/*fasta; do
	echo "Filename: $fn"
	n_lines=$(tail $fn | wc -m)

	echo "Number of lines: $n_lines"
	if [[ $n_lines -gt $threshold ]]; then
		cat $fn >> $out_dir/g3_input.fasta
	fi
done

