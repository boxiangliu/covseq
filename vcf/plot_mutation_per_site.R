library(data.table)
library(ggplot2)
library(cowplot)

in_fn = "../processed_data/vcf/count_mutation_per_site/variant_count.tsv"
out_dir = "../processed_data/vcf/plot_mutation_per_site/"
if (!dir.exists(out_dir)) {dir.create(out_dir)}

variant_count = fread(in_fn)
p = ggplot(variant_count[variant_count > 1], aes(x=pos)) + 
	geom_histogram(bins=100) + 
	xlab("SARS-CoV-2 Genomic Coordinate") + 
	ylab("Frequency of Multi-Allelic Sites") +
    theme_classic()
save_plot(sprintf("%s/distribution_of_multiallelic_sites.pdf", out_dir), p, base_height=4, base_width=6)