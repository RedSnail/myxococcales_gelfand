#!/usr/bin/env Rscript
library(ggplot2)
fam_dirs <- list.files(pattern = "*_pan$")

for (fam_dir in fam_dirs) {
  sum_file <- list.files(fam_dir, pattern="summary.txt")[[1]]
  data <- read.csv(paste0(fam_dir, "/", sum_file))
  data$sum_quali <- factor(data$sum_quali, levels=1:max(data$sum_quali))
  ggplot(data) + aes(x=sum_quali) + geom_bar(fill = "red", color="black")
  ggsave(paste(fam_dir, "hist.png", sep = "/"))
}

