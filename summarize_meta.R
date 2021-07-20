#!/usr/bin/env Rscript

genomes_table <- read.table("myxococcales_info.tsv", header=T, sep="\t")
colnames(genomes_table) <- c("assembly_level", "AC", "len", "organism", "taxid")
genomes_table$len <- round(genomes_table$len/1000000, digits=1)
fam_table <- read.table("families.csv", header=F, sep=" ", col.names=c("family", "fam_id"))

table <- cbind(genomes_table, fam_table)

write.table(table, "myxococcales_table.tsv", sep="\t", col.names=T, row.names=F, quote=F)

library(dplyr)

table %>%
  group_by(family, fam_id) %>%
  summarise(n=n(), med=median(len), min=min(len), max=max(len)) %>%
  as.data.frame() -> order_summary

order_summary <- order_summary[order(order_summary$n, decreasing=T),]
write.table(order_summary, "myxococcales_summary.tsv", sep="\t", col.names=T, row.names=F, quote=F)

