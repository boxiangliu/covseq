library(ggplot2)
library(gggenes)
library(data.table)
setwd("~/svail/open_reading_frame/scripts")

pred = fread("../processed_data/glimmer/COVID19.iterated.predict", skip=1)
setnames(pred, c("gene", "start", "end", "frame", "confidence"))
pred$molecule = "pred"
pred$confidence = NULL
pred = pred[pred$frame >= 0]

known = fread("../data/COVID19.nh")
setnames(known, c("gene", "start", "end", "frame", "name"))
known$molecule = "known"
known$name = NULL

pred[, id:=paste(start, end, sep=":")]
known[, id:=paste(start, end, sep=":")]
pred[, overlap:=id %in% known$id]
known[, overlap:=id %in% pred$id]
all = rbind(pred, known)

ggplot(all, aes(xmin=start, xmax=end, y=molecule, fill=overlap, label=gene)) + 
  geom_gene_arrow() + 
  geom_gene_label() + 
  facet_wrap(~molecule, scales="free", ncol=1)+
  scale_fill_brewer(palette="Set3")

