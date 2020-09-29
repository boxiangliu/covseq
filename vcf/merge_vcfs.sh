vcf_dir=$1
out_dir=$2
ref_fn=$3

echo "###################"
echo "# Merge VCF files #"
echo "###################"
echo "VCF: $vcf_dir"
echo "Output: $out_dir"
echo "Reference: $ref_fn"

mkdir -p $out_dir


# Remove old VCFs:
if [[ -f $out_dir/merged.vcf.gz ]]; then 
	echo "Removing outdated VCFs."
	rm $out_dir/merged.vcf.gz*
fi


# Merge VCFs:
echo "Merging VCFs."
split -a 5 -l 200 $vcf_dir/filtered_vcf.txt subset_vcfs
for i in subset_vcfs*; do
	bcftools merge -0 -l $i -Oz -o $out_dir/merged.$i.1.vcf.gz
	tabix -p vcf $out_dir/merged.$i.1.vcf.gz
done
rm subset_vcfs*

ls $out_dir/merged.*.vcf.gz > merged_vcfs
split -a 5 -l 200 merged_vcfs subset_vcfs
for i in subset_vcfs*; do
    bcftools merge -0 -l $i -Oz -o $out_dir/merged.$i.2.vcf.gz
    tabix -p vcf $out_dir/merged.$i.2.vcf.gz
done
rm subset_vcfs*
rm $out_dir/*.subset_vcfs*.1.vcf.gz*

ls $out_dir/merged.*.2.vcf.gz > merged_vcfs
bcftools merge -0 -l merged_vcfs -Oz -o $out_dir/merged.vcf.gz
tabix -p vcf $out_dir/merged.vcf.gz
rm merged_vcfs
rm $out_dir/*.subset_vcfs*.2.vcf.gz*

# Normalize merged vcf:
echo "Normalizing merged VCF."
bcftools norm -f $ref_fn $out_dir/merged.vcf.gz -m +any -Oz -o $out_dir/merged.norm.vcf.gz
rm $out_dir/merged.vcf.gz
mv $out_dir/merged.norm.vcf.gz $out_dir/merged.vcf.gz
tabix -p vcf $out_dir/merged.vcf.gz