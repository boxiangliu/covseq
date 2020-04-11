vcf_dir=../processed_data/fasta2vcf/filter_samples/
out_dir=../processed_data/merge_vcfs/merge_vcfs/
ref_fn=/Users/boxiang/Documents/work/Baidu/projects/covid_genome/data/reference/NC_045512.2.fasta
mkdir -p $out_dir

# Normalize individual vcf:
for fn in $vcf_dir/*.vcf; do
	echo $fn;
	bcftools norm -f $ref_fn $fn -Oz -o $fn.gz;
	tabix -p vcf $fn.gz;
done


# Merge VCFs:
ls $vcf_dir/*vcf.gz | split -l 200 - subset_vcfs
for i in subset_vcfs*; do
	bcftools merge -0 -l $i -Oz -o $out_dir/gisaid.$i.vcf.gz
	tabix -p vcf $out_dir/gisaid.$i.vcf.gz
done
rm subset_vcfs*

ls $out_dir/gisaid.*.vcf.gz > subset_vcfs
bcftools merge -0 -l subset_vcfs -Oz -o $out_dir/gisaid.vcf.gz
tabix -p vcf $out_dir/gisaid.vcf.gz
rm subset_vcfs
rm $out_dir/*.subset_vcfs*.vcf.gz*


# Normalize merged vcf:
bcftools norm -f $ref_fn $out_dir/gisaid.vcf.gz -m +any -Oz -o $out_dir/gisaid.norm.vcf.gz
rm $out_dir/gisaid.vcf.gz
mv $out_dir/gisaid.norm.vcf.gz $out_dir/gisaid.vcf.gz
tabix -p vcf $out_dir/gisaid.vcf.gz