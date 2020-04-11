in_dir=../processed_data/merge_vcfs/merge_vcfs/
out_dir=../processed_data/merge_vcfs/filter_sites/
mkdir -p $out_dir
bcftools view -m 2 -M 2 -r 1:1-29870 $in_dir/gisaid.vcf.gz -Oz -o $out_dir/gisaid.interim.vcf.gz 
gunzip -c $out_dir/gisaid.interim.vcf.gz | python merge_vcfs/filter_polya.py > $out_dir/gisaid.vcf
bgzip $out_dir/gisaid.vcf
tabix -p vcf $out_dir/gisaid.vcf.gz 
rm $out_dir/gisaid.interim.vcf.gz