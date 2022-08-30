#!/usr/bin/env Rscript
### INFO: Process estimated conversion rates
### DATE: 20.08.2022
### AUTHOR: Artem Baranovskii



p.in <- "/local/artem/Projects/Stability/Results/Slamseq/PCN.dnCaf1_Slamseq/conversions/"


#{SETUP}
# set variables
# ------------------------------ #
## Global
### Process arguments
if (!require("argparse")) install.packages("argparse") else library(tidyverse)
suppressPackageStartupMessages(library("argparse"))
parser$add_argument("-v", "--verbose", action = "store_true", default = TRUE,
                    help = "Print extra output [default]")






## First grab the name and extract a list of stems - for this I need one variable - in_path
p.in <- "/local/artem/Projects/Stability/Results/Slamseq/PCN.dnCaf1_Slamseq/conversions"
v.bases <- c("A", "C", "T", "G")
p.in_N <- purrr::map(v.bases,
                     ~ list.files(file.path(p.in, "genomeN", .x), pattern = "*perGene*", full.names = T)) %>% 
  purrr::reduce(., c)
p.in_NN <- purrr::map(v.bases,
                      ~ list.files(file.path(p.in, "genomeNN", .x), pattern = "*perGene*", full.names = T)) %>% 
  purrr::reduce(., c)
p.in_SNP <- list.files(file.path(p.in, "SNPs"), pattern = "*perGene*", full.names = T)


## assemble path table to ease the future sample selections
list(
  purrr::map_dfr(p.in_N,
                 ~ tibble(stem = str_split_fixed(basename(.x), "perGene.", 2)[, 2] %>% str_remove(., ".txt"),
                          base = str_split_fixed(basename(.x), ".perGene", 2)[, 1] %>% str_remove(., "counts"),
                          path_base = .x)),
  purrr::map_dfr(p.in_NN,
                 ~ tibble(stem = str_split_fixed(basename(.x), "perGene.", 2)[, 2] %>% str_remove(., ".txt"),
                          conv = str_split_fixed(basename(.x), ".perGene", 2)[, 1] %>% str_remove(., "genome|counts"),
                          path_conv = .x,
                          base = str_sub(conv, 1, 1))),
  purrr::map_dfr(p.in_SNP,
                 ~ tibble(stem = str_split_fixed(basename(.x), "perGene.", 2)[, 2] %>% str_remove(., ".txt"),
                      snp = str_split_fixed(basename(.x), ".perGene", 2)[, 1] %>% str_remove(., "genome|counts"),
                      path_snp = .x,
                      base = str_sub(snp, 1, 1)))
  ) %>% 
  purrr::reduce(., left_join, by = c("stem", "base")) %>% 
  filter(is.na(snp) | (snp == conv)) -> t.p_comb


# TODO
# 1. Function to generate CRs for all possible substitutions at timepoint Zero
comp_T0NN <- function(tab_comb, base_cutoff = 14, conv_cutoff = 5) {
  tab_comb <- tab_comb %>% filter(str_detect(stem, "T0"))
  purrr::pmap_dfr(list(tab_comb$stem,
                       tab_comb$path_base,
                       tab_comb$path_conv,
                       tab_comb$conv), 
                  function(x, y, z, r) {
                    sample_name = str_split_fixed(x, "_trimmed", 2) %>% `[`(, 1)
                    grp = str_split_fixed(sample_name, "_", 2) %>% `[`(, 1)
                    tpoint = str_split_fixed(sample_name, "_", 3) %>% `[`(, 2)
                    br = str_split_fixed(sample_name, "_", 3) %>% `[`(, 3)
                    t.base <- data.table::fread(y, data.table = F) %>% filter(V2 > base_cutoff)
                    colnames(t.base) <- c("gene", "base_count")
                    t.conv <- data.table::fread(z, data.table = F) %>% filter(V2 > conv_cutoff)
                    colnames(t.conv) <- c("gene", "conv_count")
                    t.merged <- inner_join(t.base, t.conv, by = "gene") %>% 
                      mutate(conv = r,
                             sample_name = sample_name, 
                             group = grp, 
                             timepoint = tpoint,
                             b.rep = br) %>% 
                      dplyr::select(sample_name, group, timepoint, b.rep, conv, everything())
                    }
                  ) -> t_out
  return(t_out)
}


t.cr_T0NN <- comp_T0NN(t.p_comb)



## test how it looks like?
ggplot(data = t.cr_T0NN %>% 
         mutate(col_var = ifelse(conv == "TC", "red", "black")) %>% 
         filter((conv_count / base_count) *100 < 6), 
       aes(x = conv, y = (conv_count / base_count) *100, color = col_var)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_color_manual(values = c("grey30", "red"), guide = "none") +
  stat_summary(fun.data = function(x) tibble(y = -0.85, 
                                             label = round(mean(x), 2)), 
               geom = "text", 
               vjust = 1, 
               hjust = -0.0025, 
               angle = 90) +
  geom_hline(yintercept = 0) +
  #lims(y = c(0, 10)) +
  facet_grid(group ~ b.rep) +
  #ggforce::facet_zoom(ylim = c(0, 0.6), zoom.size = , zoom.data = ifelse(conversion != "T->C", NA, FALSE), horizontal = FALSE) +
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   colour = c(rep("black", 10), "red", "black"))) +
  labs(y = "Conversion rate (%)", 
       x = "") -> check_plot


## TO zero table
### Load paths to test tables
test.N <- p.in_N[[1]]
test.NN <- p.in_NN[[1]]
test.SNP <- p.in_SNP[[1]]

### Load test tables
t.test.N <- data.table::fread(test.N)
t.test.NN <- data.table::fread(test.NN)


### figure out how to sort overlapping genes



# 2. Function to generate T -> C CRs for all timepoints


