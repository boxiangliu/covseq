library(data.table)

in_fn = "../processed_data/vcf/filter_samples/mutation_count.tsv"
mutation_count = fread(in_fn, sep="\t", col.names=c("sample", "mutation"))
p = ggplot(mutation_count, aes(x=mutation)) + geom_histogram(bins=300) + scale_y_sqrt() + scale_x_log10(breaks=c(1,10,100,150,1000,10000)) + geom_vline(xintercept=150, color="red", linetype="dashed")
print(p)