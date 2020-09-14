library(data.table)
library(ggplot2)
library(cowplot)
in_fn = "../processed_data/preprocess/filter_fasta/genome_lengths.tsv"
out_dir = "../processed_data/preprocess/plot_genome_lengths/"
if (!dir.exists(out_dir)) {dir.create(out_dir)}

genome_lengths = fread(in_fn, sep="\t", col.names = c("description", "length"))
p = ggplot(genome_lengths, aes(x=length)) + 
	geom_histogram(bins=100) + 
	scale_y_sqrt() + 
	geom_vline(xintercept = 25000, color="red", linetype="dashed") + 
	scale_x_continuous(breaks=c(0,5000,10000,15000,20000,25000,30000)) +
	xlab("SARS-CoV-2 Sequence Lengths") + 
	ylab("Frequency") + 
    theme_classic()
print(p)

save_plot(sprintf("%s/genome_lengths.pdf", out_dir), p, base_width=6, base_height=3)