fasta_dir=../data/GISAID_Human_Coronavirus_fasta_clean/
vapid_dir=../src/vapid/
out_dir=../processed_data/parse_vapid/
mkdir -p $out_dir/{vapid,parsed}

cat $fasta_dir/EPI_ISL_*.fasta > $fasta_dir/all.fasta

echo "strain,collection-date,country,coverage" > $fasta_dir/metadata.csv
for f in $fasta_dir/EPI_ISL_*.fasta; do
	stem=$(basename $f .fasta)
	echo "${stem},2020,World,0" >> $fasta_dir/metadata.csv
done

module load mafft
python3 $vapid_dir/vapid3.py \
	--db $vapid_dir/all_virus.fasta \
	--r NC_045512.2 \
	--metadata_loc $fasta_dir/metadata.csv \
	$fasta_dir/all.fasta \
	$vapid_dir/template.sbt
mv EPI_ISL_* $out_dir/vapid/

for strain in $out_dir/vapid/EPI_ISL_*/; do
	base=$(basename $strain)
	# echo $strain/$base
	python parse_vapid.py --in_prefix $strain/$base --out_dir $out_dir/parsed/
done


