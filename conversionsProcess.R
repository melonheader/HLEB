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
# 1. Function to generate CRs for all possible substitutions at timepoint Zero ---- DONE
# 1.1. Function to plot CRs at timepoint zero ---- DONE
# 2. Function to read in all T -> C conversions, don't forget cutoffs ---- DONE
# 2.1.Function to filter T -> C conversions ---- DONE
# 2.2. Function to plot CRs as boxes per TP per replicate ---- DONE
# 2.3. Function to detect otlying replicates (KS test) ---- DONE
# 3. Function to normalize convertions to timepoint zero and compute averages ---- DONE
# 4. Write a script to estimate Half-lives ---- MORE EFFICIENT ---- break into bulks of trnascripts and parallelize



















## workflow
t.cr_TC <- comp_TC(t.p_comb)
plot_TC.qc_boxes(t.cr_TC)
plot_TC.qc_ks(t.cr_TC)
# Drop 
# GFP: T0 br1, 
# dnCaf1: T0 br3, T16 br 2
filt_out <- filt_TC(t.cr_TC, 
                    min_T0 = 2, 
                    min_T4 = 2, 
                    min_T8 = 0,
                    min_T16 = 0,
                    out_drop = c("GFP_T0_br1", "dnCaf1_T0_br1", "dnCaf1_T16_br2"))
t.cr_TC.filt <- filt_out$tab_cr_TC.filt
t.cr_TC.filt.norm <- norm_TC(t.cr_TC.filt)$normalised









## Try two strategies ----
#### first normalize, then average
### test how many genes we will have in approach number one ---- NE IMEET SMISLAs
t.cr_TC.filt %>% 
  dplyr::select(-sample_name, -contains("count"), -cr_perc, -conv) %>% 
  mutate(n = 1) %>% 
  pivot_wider(., names_from = timepoint, values_from = n, values_fill = 0) %>% 
  mutate(n_sum = T0 + T4 + T8 + T16) -> rep_test
rep_test %>% filter(T0 != 0) %>%  tabyl(b.rep, n_sum, group) ### looks surprisingly okay ---- proceed with 3 TP minimum



## Overall, USE A NOTEBOOK











