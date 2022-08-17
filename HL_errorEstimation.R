
######################################
### DATE: 15.07.2022
### INFO: Estimate error of half-lives
######################################



# {SETUP}
## Paths
p.parent.dir <- getwd()

if (!require("tidyverse")) install.packages("tidyverse") else library(tidyverse)
if (!require("magrittr")) install.packages("magrittr") else library(magrittr)
if (!require("purrr")) install.packages("purrr") else library(purrr)
if (!require("tictoc")) install.packages("tictoc") else library(tictoc)






source("~/Scripts/Notebooks/Stability/Reviews/F_HL_errorEstimation.R")


## here goes genomic conversion tables processing....






## First write a template script using compartment separation slamseq half-lives
d.s <- read_tsv("/local/artem/Projects/Stability/Results/hl_estimation/Slamseq/data_s.tsv") %>% 
  dplyr::rename(Group = Compartment)

estimate_hl(project = "Stability",
            experiment = "Compartments_HlErrors",
            data_s = d.s,
            time_limit = 24, 
            warnings = F,
            aux_results = T) -> t.hl_comp.Err
write_tsv(t.hl_comp.Err, "/local/artem/Projects/Stability/Results/hl_estimation/SlamseqComp/data_hl.tsv")

















## Process table from slamdunk and estimate half-lives
# l.cvs <- list.files("/local/artem/Projects/Stability/Results/PCN.Slamseq_dnCaf1/slamdunk/count", 
#                     pattern = "*_collapsed.csv", 
#                     full.names = T)
# t.cvs <- purrr::map(l.cvs,
#                     ~ data.table::fread(.x) %>% 
#                       as_tibble() %>% 
#                       mutate(Sample = str_split_fixed(basename(.x), "_trimmed", 2)[, 1])
#                     ) %>% 
#   purrr::reduce(., rbind)
# 
# ## Filter
# ## readsCPM
# 
# # global parameters
# T_filter <- 10
# TC_filter <- 0
# 
# 
# # Sample, Group, Timepoint, Replicate, Name, T_Count, TC_Count
# prepare_sdunk_tcount <- function(t.tcount, 
#                                  r_filter = 2, 
#                                  dr_filter = 2, 
#                                  T_filter = NULL, 
#                                  TC_filter = NULL, 
#                                  stddev_norm_cutoff = NULL, 
#                                  mm10_gene_names = NULL) {
#   ## 
#   t.tcount <- t.tcount %>% 
#     mutate(Group = str_split_fixed(Sample, "_", 2)[, 1],
#            Timepoint = str_split_fixed(Sample, "_", 3)[, 2],
#            Replicate = str_split_fixed(Sample, "_", 3)[, 3]) %>% 
#     dplyr::select(Sample, Group, Timepoint, Replicate, Name = gene_name, length, readsCPM, T_count = coverageOnTs, TC_count = conversionsOnTs, Conv_Rate = conversionRate)  
#     
#   # T_count filter
#   message("filtering by T_count...")
#   if (purrr::is_null(T_filter)) {
#     stop("Specify T_filter")
#   }
#   data_filtered <- t.tcount %>% 
#     filter(T_count > T_filter)
#   
#   # TC_count filter
#   message("filtering by TC_count...")
#   if (purrr::is_null(TC_filter)) {
#     stop("Specify TC_filter")
#   }
#   tr_to_drop <- data_filtered %>% 
#     filter(str_detect(Timepoint, "0"), 
#            TC_count < TC_filter) %>% 
#     distinct(Name) %>% 
#     pull(Name)
#   data_filtered <- filter(data_filtered, !Name %in% tr_to_drop)
#   
#   #
#   # filter transcripts undetected at T0 - mandatory
#   message("filtering timepoint 0 dropouts...")
#   r_check <- data_filtered %>%
#     dplyr::select(Name, Sample, Conv_Rate) %>%
#     myspread(., 
#              key = Sample, 
#              value = Conv_Rate) %>%
#     gather(., 
#            contains("Conv_Rate"), 
#            key = "Sample", 
#            value = "Conv_Rate") %>%
#     mutate(Group = str_split_fixed(Sample, "_", Inf)[, 1], 
#            Timepoint = str_split_fixed(Sample, "_", Inf)[, 2], 
#            Replicate = str_split_fixed(Sample, "_", Inf)[, 3]) %>%
#     group_by(Name, Group, Timepoint) %>% 
#     summarise(obs_count = sum(!is.na(Conv_Rate)), 
#               .groups = "drop_last")
#   tr_to_drop <- r_check %>% 
#     filter(str_detect(Timepoint, "0"), 
#            obs_count < r_filter) %>% 
#     distinct(Name) %>% 
#     pull(Name)
#   data_filtered <- filter(data_filtered, !Name %in% tr_to_drop)
#   
#   
#   # Omit transripts where more than 2 timepoints went undetected - mandatory
#   message("filtering dropouts across all timepoints...")
#   dr_check <- list()
#   for (i in names(table(data_filtered$Group))) {
#     dr_table <- tibble()
#     dat_tmp <- filter(data_filtered, Group == i)
#     for (y in unique(data_filtered$Name)) {
#       obs_count <- dat_tmp %>% 
#         filter(Name == y) %>% 
#         group_by(Group, Timepoint, Name) %>% 
#         summarise(count = n(), 
#                   .groups = "drop_last") %>% 
#         dim(.) %>%
#         `[`(1)
#       dr_tmp <- tibble(Name = y,
#                        Group = i, 
#                        Dr_count = 4 - obs_count)
#       dr_table <- rbind(dr_table, dr_tmp)
#     }
#     dr_check[[i]] <- dr_table
#   }
#   dr_check <- purrr::reduce(dr_check, rbind)
#   tr_to_drop <- dr_check %>% 
#     filter(Dr_count >= dr_filter) %>% 
#     distinct(Name) %>% 
#     pull(Name)
#   data_filtered <- filter(data_filtered, !Name %in% tr_to_drop)
#   
#   # load gene names
#   if (purrr::is_null(mm10_gene_names)) {
#     mm10_gene_names <- read_tsv("/local/artem/Data/Annotations/mm10_gene_names.tsv", 
#                                 col_names = TRUE,
#                                 col_types = cols(.default = "c"))
#   }
#   
#   # prepare table
#   data_s <- data_filtered %>% 
#     inner_join(mm10_gene_names, ., 
#                by = c("gene_id" = "Name")) %>% 
#     mutate(Conv_Rate = TC_count / T_count, 
#            Timepoint = as.numeric(as.character(str_remove_all(Timepoint, "T|hrs"))))
#   
#   # estimate and normalize conversion rates
#   ## compute errors and means for t0
#   message("estimating means and errors for timepoint 0...")
#   data_t0 <- list()
#   for (i in names(table(data_s$Group))) {
#     ### define
#     t0 <- tibble()
#     ### subset group
#     t0_table <- data_s %>% 
#       filter(Group == i & Timepoint == 0)
#     for (y in unique(data_s$gene_id)) {
#       ### subset gene
#       t0_sample <- t0_table %>% 
#         filter(gene_id == y) %>% 
#         pull(Conv_Rate)
#       t0_stddev <- sd(t0_sample)
#       t0_stddev_norm <- t0_stddev / mean(t0_sample)
#       t0_tmp <- tibble(gene_id = y, 
#                        Group = i, 
#                        cr_mean = mean(t0_sample), 
#                        stddev = t0_stddev,
#                        stddev_norm = t0_stddev_norm)
#       t0 <- rbind(t0, t0_tmp)
#     }
#     data_t0[[i]] <- t0
#   }
#   
#   ## apply cutoff based on t0 error
#   if (!purrr::is_null(stddev_norm_cutoff)) {
#     g_passed <- list()
#     for (i in names(table(data_s$Group))) {
#       g_passed[[i]] <- data_t0[[i]] %>% 
#         filter(stddev_norm < stddev_norm_cutoff) %>% 
#         pull(gene_id)
#     }
#     g_passed <- purrr::reduce(g_passed, base::intersect)
#     data_s <- data_s %>% 
#       filter(gene_id %in% g_passed)
#     data_t0 <- map(data_t0, ~ filter(., gene_id %in% g_passed))
#   }
#   
#   ## compute normalized conversion rates
#   message("estimating normalized conversion rates...")
#   data_s$Conv_Rate_Norm <- 0
#   data_s <- data_s %>% 
#     mutate(Conv_Rate_Norm = ifelse(Timepoint == 0, 1, Conv_Rate_Norm))
#   for (i in names(table(data_s$Group))) {
#     for (y in unique(data_s$gene_id)) {
#       data_s <- data_s %>% 
#         mutate(Conv_Rate_Norm = ifelse(Timepoint != 0 & Group == i & gene_id == y, 
#                                        round(Conv_Rate / data_t0[[i]][data_t0[[i]]$gene_id == y, ][["cr_mean"]], 3), 
#                                        Conv_Rate_Norm
#         )
#         )
#     }
#   }
#   
#   ## artificially set norm_conversion rates above 1 to 1 - we cannot adequately estimate the error 
#   data_s <- data_s %>% 
#     mutate(Conv_Rate_Norm = ifelse(Conv_Rate_Norm > 1, 1, Conv_Rate_Norm))
#   
#   ## compute means & errors for other timepoints 
#   message("estimating means and errors for normalized conversion rates...")
#   data_s_s <- data_s %>% 
#     group_by(gene_id, gene_name, Timepoint, Group) %>% 
#     summarise(Mean_Conv_Rate_Norm = mean(Conv_Rate_Norm), 
#               Stddev_Norm = sd(Conv_Rate_Norm), 
#               .groups = "drop_last") %>% 
#     ungroup()
#   for (i in names(table(data_s_s$Group))) {
#     for (y in unique(data_s_s$gene_id)) {
#       data_s_s <- data_s_s %>% 
#         mutate(Stddev_Norm = ifelse(Timepoint == 0 & Group == i & gene_id == y, 
#                                     data_t0[[i]][data_t0[[i]]$gene_id == y, ][["stddev_norm"]], 
#                                     Stddev_Norm)
#         )
#     }
#   }
#   # Enlist
#   output <- list(data_s_s, data_s)
#   names(output) <- c("data_s_s", "data_s")
# }



