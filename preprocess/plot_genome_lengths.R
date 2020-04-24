in_fn = "../../processed_data/preprocess/filter_fasta/genome_lengths.tsv"
genome_lengths = fread(in_fn, sep="\t", col.names = c("description", "length"))
p = ggplot(genome_lengths, aes(x=length)) + geom_histogram(bins=100) + scale_y_sqrt() + geom_vline(xintercept = 25000, color="red", linetype="dashed") + geom_text(x = 25000, y = 30, label = "Cutoff = 25000", color="red", hjust=1.1, size = 5)
print(p)