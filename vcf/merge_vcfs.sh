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
split -l 200 $vcf_dir/filtered_vcf.txt subset_vcfs
for i in subset_vcfs*; do
	bcftools merge -0 -l $i -Oz -o $out_dir/merged.$i.vcf.gz
	tabix -p vcf $out_dir/merged.$i.vcf.gz
done
rm subset_vcfs*

ls $out_dir/merged.*.vcf.gz > subset_vcfs
bcftools merge -0 -l subset_vcfs -Oz -o $out_dir/merged.vcf.gz
tabix -p vcf $out_dir/merged.vcf.gz
rm subset_vcfs
rm $out_dir/*.subset_vcfs*.vcf.gz*


# Normalize merged vcf:
echo "Normalizing merged VCF."
bcftools norm -f $ref_fn $out_dir/merged.vcf.gz -m +any -Oz -o $out_dir/merged.norm.vcf.gz
rm $out_dir/merged.vcf.gz
mv $out_dir/merged.norm.vcf.gz $out_dir/merged.vcf.gz
tabix -p vcf $out_dir/merged.vcf.gz