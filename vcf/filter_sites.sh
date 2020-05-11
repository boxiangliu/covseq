in_fn=$1
out_prefix=$2

# in_dir=../processed_data/merge_vcfs/merge_vcfs/
# out_dir=../processed_data/merge_vcfs/filter_sites/

echo "######################"
echo "# Filter merged VCF  #"
echo "######################"
echo "Input VCF: $in_fn"
echo "Output VCF prefix: $out_prefix"

if [[ -f $out_prefix.vcf.gz ]]; then
	echo "Removing old VCFs."
	rm $out_prefix.vcf.gz{,.tbi}
fi

bcftools view -m 2 -M 2 $in_fn -Ov -o interim.vcf 
cat interim.vcf | python vcf/filter_polya.py > $out_prefix.vcf
rm interim.vcf
bgzip $out_prefix.vcf
tabix -p vcf $out_prefix.vcf.gz
