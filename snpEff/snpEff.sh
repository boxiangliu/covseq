in_fn=$1
out_prefix=$2
echo "##################"
echo "# Annotating VCF #"
echo "##################"
echo "Input: $in_fn"
echo "Output: $out_prefix"

java -jar ext/snpEff/snpEff.jar NC_045512.2 $in_fn > $out_prefix.vcf
bgzip $out_prefix.vcf
tabix -p vcf $out_prefix.vcf.gz
rm snpEff_genes.txt snpEff_summary.html
