#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
dir_prefix <- args[[1]]
alias <- substr(dir_prefix, 1, 4)

venn <- "venn" %in% args
plotlist <- list()

library(tidyverse)
library(gsubfn)

metadata <- read_tsv("myxococcales_table.tsv")

pantable_file <- fn$Sys.glob("`dir_prefix`_pan/PanGenome-`alias`-*.lst.quanti.txt")[[1]]
pantable <- fn$read_csv(pantable_file)
pantable %>%
  mutate(name=fam_num) %>%
  select(-fam_num) -> pantable

is_zero <- function(x) {
  sum(x) > 0
}

library(ggVennDiagram)

if (venn) {
  big_mainstream <- c("Chondromyces", "Sorangium", "Minicystis", "Labilithrix")
  big_altrock <- c("Corallococcus", "Myxococcus", "Stigmatella", "Melittangium", "Cystobacter", "Archangium")
  
  small_mainstream <- c("Sandaracinus", "Haliangium")
  small_altrock <- c("Vulgatibacter", "Anaeromyxobacter")
  
  clade_list <- list()
  clade_list[c(big_mainstream, small_mainstream)] <- "mainstream"
  clade_list[c(big_altrock, small_altrock)] <- "altrock"
  
  
  size_list <- list()
  size_list[c(big_altrock, big_mainstream)] <- "big"
  size_list[c(small_mainstream, small_altrock)] <- "small"
  
  metadata %>%
    select(organism, AC, len) %>%
    mutate(genus=word(organism, 1)) %>%
    filter(genus %in% names(clade_list)) %>%
    mutate(size=unlist(size_list[genus]), clade=unlist(clade_list[genus])) %>%
    left_join(fn$read_tsv("`dir_prefix`_pan/names_mapping.tsv", col_names = c("AC", "name", "organism")), by=c("AC", "organism")) -> meta_table
  
  meta_table %>%
    left_join(pantable, by="name") -> common_table
  
  common_table %>%
    group_by(clade, size) %>%
    summarise_at(vars(matches("[0-9]")), is_zero) -> summary_result
  
  summary_result %>%
    pivot_longer(cols = matches("[0-9]"), names_to="family") %>%
    filter(value) %>%
    group_by(clade, size) %>%
    select(family, clade, size) %>%
    group_split(.keep = F) %>%
    lapply(function(x) {
      x$family
    }) -> set_list
  
  venn_data <- process_data(Venn(set_list))
  
  ggplot() +
    geom_sf(aes(fill=count), data = venn_region(venn_data)) +
    geom_sf(aes(color=name), data = venn_setedge(venn_data), show.legend = F) +
    geom_sf_text(aes(label = name), data = venn_setlabel(venn_data) %>% 
                   mutate(name=c("big altrock", "small altrock", "big mainstream", "small maistream"))) +
    geom_sf_text(aes(label=count), fontface = "bold", data = venn_region(venn_data)) + 
    scale_fill_gradient(low = "white", high = "red") + 
    scale_color_discrete(type = c("darkgreen", "lightgreen", "darkblue", "lightblue")) -> venn_plot
  
  plotlist$Venn_plots <- venn_plot
}

pantable %>%
  pivot_longer(cols = matches("[0-9]"), names_to="family") -> longer_pantable

longer_pantable %>%
  group_by(value) %>%
  summarise(quanti=n()) -> copyn_distribution

ggplot(copyn_distribution, aes(x=value, y=quanti)) + 
  geom_bar(stat = "identity") +
  scale_y_log10() +
  ggtitle("Histogram of quantity of members of familily in genome") +
  xlab("members") + 
  ylab("quantity") -> copyn_hist

plotlist$Copyn_hists <- copyn_hist



longer_pantable %>%
  filter(value>0) %>%
  group_by(family) %>%
  summarise(genomes=n()) %>%
  group_by(genomes) %>%
  summarise(families=n()) -> family_by_genomes_shared_distr

ggplot(family_by_genomes_shared_distr, aes(x=genomes, y=families)) +
  geom_bar(stat = "identity") +
  ggtitle("Distribution of families by genomes which share it") +
  xlab("quantity of genomes") +
  ylab("quantity of families") -> fam_distr_plot

plotlist$Fam_distrs <- fam_distr_plot


# fair evgeniy gordienko plots

gen_fam_distr <- function(data, replics) {
  data %>%
    mutate(value = value>0) %>%
    group_by(family) %>%
    summarise(presence = sum(value)/n()) %>%
    filter(presence>0) %>%
    mutate(presence=cut(presence, 
                        breaks = c(0, 0.1, 0.5, 0.9, 1, Inf), 
                        right = F, 
                        labels = c("sparse cloud", "mid-cloud", "dense cloud", "soft core", "core"))) %>%
    group_by(presence, .drop=F) %>%
    summarise(nfam = n())  
}

nrepls = 10
repl_ngenomes_sampling <- function(ngenomes) {
  bind_rows(replicate(nrepls, 
                      longer_pantable %>%
                        group_by(name) %>%
                        nest() %>%
                        ungroup() %>%
                        sample_n(ngenomes) %>%
                        unnest(cols = c(data)), 
                      simplify = F), 
            .id = "replics") %>%
    group_by(replics) %>%
    group_modify(gen_fam_distr) %>%
    group_by(presence) %>%
    summarise(mean_fams = mean(nfam), sdev = sd(nfam)) %>%
    mutate(ngenomes=ngenomes)
  
}

bind_rows(lapply(1:nrow(pantable), repl_ngenomes_sampling)) -> permuts

library(spatstat.utils)

permuts %>%
  group_by(ngenomes) %>%
  group_modify(function(.x, .y) {
    .x %>%
      mutate(ycumsum=revcumsum(mean_fams))
  }) -> for_error_plot

ggplot(for_error_plot, aes(x=ngenomes, y=mean_fams, color=presence)) +
  geom_point(stat = "identity", position = "stack") +
  geom_line(position = "stack") +
  geom_errorbar(aes(ymin=ycumsum-sdev, ymax=ycumsum+sdev)) +
  ggtitle(paste0("Eugene Gordienko plot with error bars. npermuts: ", nrepls)) +
  xlab("genomes taken") +
  ylab("quantity of genes") -> error_plot

plotlist$Error_plots <- error_plot



# rigged evgeniy gordienko plots, with expectetions instead of fair permutations
longer_pantable %>%
  mutate(value = value>0) %>%
  group_by(family) %>%
  summarise(presence = sum(value)/n()) -> family_fractions

bind_rows(lapply(1:nrow(pantable), function(n) {
  bind_rows(lapply(n:1, function(j) {
    family_fractions %>%
      mutate(expectation=choose(n, j)*(presence^j)*((1 - presence)^(n-j))) %>%
      summarise(sum_exp=sum(expectation)) %>%
      mutate(nsample=n, fraction=j/n)
  }))
})) %>%
  mutate(occurence_group=cut(fraction, 
                      breaks = c(0, 0.1, 0.5, 0.9, 1, Inf), 
                      right = F, 
                      labels = c("sparse cloud", "mid-cloud", "dense cloud", "soft core", "core"))) -> expected_counts

expected_counts %>%
  group_by(nsample, occurence_group) %>%
  summarise(total_count=sum(sum_exp)) %>%
  group_by(nsample) %>%
  group_modify(function(.x, .y) {
    mutate(.x, cumsums=revcumsum(total_count))
  }) -> by_occ_group

ggplot(expected_counts, aes(x=nsample, y=sum_exp, fill=fraction)) +
  geom_bar(stat = "identity") + 
  scale_fill_gradient(low = "green", high = "red") +
  geom_segment(aes(color=occurence_group, y=cumsums, yend=cumsums, xend=nsample+0.5, x=nsample-0.5), data=by_occ_group, inherit.aes = F) +
  scale_color_grey() +
  xlab("genomes taken") +
  ylab("quantity of genes") + 
  ggtitle("Eugene Gordienko plot") -> gord_plot

plotlist$Gord_plots <- gord_plot

for (dir_to_save in names(plotlist)) {
  plot <- plotlist[[dir_to_save]]
  
  dir.create(dir_to_save, showWarnings = FALSE)
  picname <- paste(dir_to_save, paste(dir_prefix, "png", sep = "."), sep = "/")
  ggsave(picname, plot + labs(subtitle = dir_prefix))
}

