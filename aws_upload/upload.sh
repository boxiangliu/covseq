# Zip up VCF files and metadata files:
wd=$(pwd)
cd ../data/aggregated/vcf/merged/
rm VCF.zip
zip VCF.zip *.{gz,gz.tbi} README
cd $wd
cd ../data/aggregated/metadata/
rm metadata.zip
zip metadata.zip *.tsv README
cd $wd 

# Upload to AWS
aws s3 cp ../data/aggregated/vcf/merged/VCF.zip s3://covseq/VCF/ --acl public-read
aws s3 cp ../data/aggregated/metadata/metadata.zip s3://covseq/metadata/ --acl public-read
