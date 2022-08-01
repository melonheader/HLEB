### INFO: Half-life estimation
### DATE: 12.09.2019
### AUTHOR: Artem Baranovskii
library(tidyverse)
library(purrr)
library(foreach)

# ---------------------------------------------------------------------------- #
load_TC_counts <- function(data_path,
                           time_pts = c("0hrs", "4hrs", "8hrs", "16hrs"),
                           name_sub = c(18, -5),
                           grp_splt = c(3, 1),
                           tp_splt = c(3, 2),
                           rep_sub = c(-7, -1),
                           recount_test = TRUE) {
  ##Read in tcounts tables
  if (recount_test == TRUE) {
    #
    t_path <- file.path(data_path, "genomeN/T/")
    tc_path <- file.path(data_path, "genomeNN/T/")
    #
    t_counts_list <- list.files(path = t_path, pattern = 'countsT.perGene.', full.names = TRUE) %>%
      purrr::map(~ read_tsv(file = .,
                            col_names = FALSE,
                            col_types = cols(
                              X1 = col_character(),
                              X2 = col_double()))) %>%
      purrr::map(~ `colnames<-`(., c("Name", "T_Count"))) %>%
      purrr::map2(., str_sub(list.files(path = t_path,
                                        pattern = 'countsT.perGene.'),
                             (name_sub[1] - 1),
                             name_sub[2]),
                  ~ mutate(.x, Sample = paste0(.y))) %>%
      set_names(., str_sub(list.files(path = t_path, pattern = 'countsT.perGene.'), (name_sub[1] - 1) , name_sub[2]))
    #
    tc_counts_list <- list.files(path = tc_path, pattern = 'genomeTC.perGene.', full.names = TRUE) %>% 
      purrr::map(~ read_tsv(file = .,
                            col_names = FALSE, 
                            col_types = cols( 
                              X1 = col_character(),
                              X2 = col_double()))) %>% 
      purrr::map(~ `colnames<-`(., c("Name", "TC_Count"))) %>% 
      purrr::map2(., str_sub(list.files(path = tc_path, 
                                        pattern = 'genomeTC.perGene.'), 
                             name_sub[1], 
                             name_sub[2]), 
                  ~ mutate(.x, Sample = paste0(.y))) %>%
      set_names(., str_sub(list.files(path = tc_path, pattern = 'genomeTC.perGene.'), name_sub[1], name_sub[2]))
    
  } else {
    #
    t_counts_list <- list.files(path = file.path(data_path, "genomeT/"), pattern = 'countsT.perGene.', full.names = TRUE) %>% 
      purrr::map(.,
                 ~ read_tsv(file = .,
                            col_names = FALSE, 
                            col_types = cols(
                              X1 = col_character(),
                              X2 = col_double()))) %>% 
      purrr::map(~ `colnames<-`(., c("Name", "T_Count"))) %>% 
      purrr::map2(., str_sub(list.files(path = file.path(data_path, "genomeT/"), 
                                        pattern = 'countsT.perGene.'), 
                             (name_sub[1] - 1), 
                             name_sub[2]), 
                  ~ mutate(.x, Sample = paste0(.y))) %>%
      set_names(., str_sub(list.files(path = file.path(data_path, "genomeT/"), pattern = 'countsT.perGene.'), (name_sub[1] - 1) , name_sub[2]))
    #
    tc_counts_list <- list.files(path = file.path(data_path, "genomeTC/"), pattern = 'genomeTC.perGene.', full.names = TRUE) %>% 
      purrr::map(~ read_tsv(file = .,
                            col_names = FALSE, 
                            col_types = cols( 
                              X1 = col_character(),
                              X2 = col_double()))) %>% 
      purrr::map(~ `colnames<-`(., c("Name", "TC_Count"))) %>% 
      purrr::map2(., str_sub(list.files(path = file.path(data_path, "genomeTC/"), 
                                        pattern = 'genomeTC.perGene.'), 
                             name_sub[1], 
                             name_sub[2]), 
                  ~ mutate(.x, Sample = paste0(.y))) %>%
      set_names(., str_sub(list.files(path = file.path(data_path, "genomeTC/"), pattern = 'genomeTC.perGene.'), name_sub[1], name_sub[2]))
  }
  data_list <- inner_join(purrr::reduce(t_counts_list, rbind), 
                          purrr::reduce(tc_counts_list, rbind), 
                          by = c("Sample", "Name")) %>% 
    mutate(Group = str_split_fixed(Sample, "_", grp_splt[1])[, grp_splt[2]], 
           Timepoint = factor(str_split_fixed(Sample, "_", tp_splt[1])[, tp_splt[2]], levels = time_pts),
           Replicate = str_sub(Sample, rep_sub[1], rep_sub[2])) %>%
    dplyr::select(Sample, Group, Timepoint, Replicate, Name, T_Count, TC_Count)
  return(data_list)
}



# ---------------------------------------------------------------------------- #
prepare_data <- function(data_list, 
                         r_filter = 2, 
                         dr_filter = 2, 
                         T_filter = NULL, 
                         TC_filter = NULL, 
                         stddev_norm_cutoff = NULL, 
                         mm10_gene_names = NULL) { 
  # T_count filter
  message("filtering by T_count...")
  if (purrr::is_null(T_filter)) {
    stop("Specify T_filter")
  }
  data_filtered <- data_list %>% 
    filter(T_Count > T_filter)
  
  # TC_count filter
  message("filtering by TC_count...")
  if (purrr::is_null(TC_filter)) {
    stop("Specify TC_filter")
  }
  tr_to_drop <- data_filtered %>% 
    filter(str_detect(Timepoint, "0"), 
           TC_Count < TC_filter) %>% 
    distinct(Name) %>% 
    pull(Name)
  data_filtered <- filter(data_filtered, !Name %in% tr_to_drop)
  
  # filter transcripts undetected at T0 - mandatory
  message("filtering timepoint 0 dropouts...")
  r_check <- data_filtered %>%
    mutate(Conv_Rate = round(TC_Count / T_Count, 5), 
           Sample = str_c(Group, Timepoint, Replicate, sep = "_")) %>%
    dplyr::select(Name, Sample, Conv_Rate) %>%
    myspread(., 
             key = Sample, 
             value = Conv_Rate) %>%
    gather(., 
           contains("Conv_Rate"), 
           key = "Sample", 
           value = "Conv_Rate") %>%
    mutate(Group = str_split_fixed(Sample, "_", Inf)[, 1], 
           Timepoint = str_split_fixed(Sample, "_", Inf)[, 2], 
           Replicate = str_split_fixed(Sample, "_", Inf)[, 3]) %>%
    group_by(Name, Group, Timepoint) %>% 
    summarise(obs_count = sum(!is.na(Conv_Rate)), 
              .groups = "drop_last")
  tr_to_drop <- r_check %>% 
    filter(str_detect(Timepoint, "0"), 
           obs_count < r_filter) %>% 
    distinct(Name) %>% 
    pull(Name)
  data_filtered <- filter(data_filtered, !Name %in% tr_to_drop)
  
  # Omit transripts where more than 2 timepoints went undetected - mandatory
  message("filtering dropouts across all timepoints...")
  dr_check <- list()
  for (i in names(table(data_filtered$Group))) {
    dr_table <- tibble()
    dat_tmp <- filter(data_filtered, Group == i)
    for (y in unique(data_filtered$Name)) {
      obs_count <- dat_tmp %>% 
        filter(Name == y) %>% 
        group_by(Group, Timepoint, Name) %>% 
        summarise(count = n(), 
                  .groups = "drop_last") %>% 
        dim(.) %>%
        `[`(1)
      dr_tmp <- tibble(Name = y,
                       Group = i, 
                       Dr_count = 4 - obs_count)
      dr_table <- rbind(dr_table, dr_tmp)
    }
    dr_check[[i]] <- dr_table
  }
  dr_check <- purrr::reduce(dr_check, rbind)
  tr_to_drop <- dr_check %>% 
    filter(Dr_count >= dr_filter) %>% 
    distinct(Name) %>% 
    pull(Name)
  data_filtered <- filter(data_filtered, !Name %in% tr_to_drop)
  
  # load gene names
  if (purrr::is_null(mm10_gene_names)) {
    mm10_gene_names <- read_tsv("/local/artem/Data/Annotations/mm10_gene_names.tsv", 
                                col_names = TRUE,
                                col_types = cols(.default = "c"))
  }
  
  # prepare table
  data_s <- data_filtered %>% 
    inner_join(mm10_gene_names, ., 
               by = c("gene_id" = "Name")) %>% 
    mutate(Conv_Rate = TC_Count / T_Count, 
           Timepoint = as.numeric(as.character(str_remove_all(Timepoint, "T|hrs"))))
  
  # estimate and normalize conversion rates
  ## compute errors and means for t0
  message("estimating means and errors for timepoint 0...")
  data_t0 <- list()
  for (i in names(table(data_s$Group))) {
    ### define
    t0 <- tibble()
    ### subset group
    t0_table <- data_s %>% 
      filter(Group == i & Timepoint == 0)
    for (y in unique(data_s$gene_id)) {
      ### subset gene
      t0_sample <- t0_table %>% 
        filter(gene_id == y) %>% 
        pull(Conv_Rate)
      t0_stddev <- sd(t0_sample)
      t0_stddev_norm <- t0_stddev / mean(t0_sample)
      t0_tmp <- tibble(gene_id = y, 
                       Group = i, 
                       cr_mean = mean(t0_sample), 
                       stddev = t0_stddev,
                       stddev_norm = t0_stddev_norm)
      t0 <- rbind(t0, t0_tmp)
    }
    data_t0[[i]] <- t0
  }
  
  ## apply cutoff based on t0 error
  if (!purrr::is_null(stddev_norm_cutoff)) {
    g_passed <- list()
    for (i in names(table(data_s$Group))) {
      g_passed[[i]] <- data_t0[[i]] %>% 
        filter(stddev_norm < stddev_norm_cutoff) %>% 
        pull(gene_id)
    }
    g_passed <- purrr::reduce(g_passed, base::intersect)
    data_s <- data_s %>% 
      filter(gene_id %in% g_passed)
    data_t0 <- map(data_t0, ~ filter(., gene_id %in% g_passed))
  }
  
  ## compute normalized conversion rates
  message("estimating normalized conversion rates...")
  data_s$Conv_Rate_Norm <- 0
  data_s <- data_s %>% 
    mutate(Conv_Rate_Norm = ifelse(Timepoint == 0, 1, Conv_Rate_Norm))
  for (i in names(table(data_s$Group))) {
    for (y in unique(data_s$gene_id)) {
      data_s <- data_s %>% 
        mutate(Conv_Rate_Norm = ifelse(Timepoint != 0 & Group == i & gene_id == y, 
                                       round(Conv_Rate / data_t0[[i]][data_t0[[i]]$gene_id == y, ][["cr_mean"]], 3), 
                                       Conv_Rate_Norm
        )
        )
    }
  }
  
  ## artificially set norm_conversion rates above 1 to 1 - we cannot adequately estimate the error 
  data_s <- data_s %>% 
    mutate(Conv_Rate_Norm = ifelse(Conv_Rate_Norm > 1, 1, Conv_Rate_Norm))
  
  ## compute means & errors for other timepoints 
  message("estimating means and errors for normalized conversion rates...")
  data_s_s <- data_s %>% 
    group_by(gene_id, gene_name, Timepoint, Group) %>% 
    summarise(Mean_Conv_Rate_Norm = mean(Conv_Rate_Norm), 
              Stddev_Norm = sd(Conv_Rate_Norm), 
              .groups = "drop_last") %>% 
    ungroup()
  for (i in names(table(data_s_s$Group))) {
    for (y in unique(data_s_s$gene_id)) {
      data_s_s <- data_s_s %>% 
        mutate(Stddev_Norm = ifelse(Timepoint == 0 & Group == i & gene_id == y, 
                                    data_t0[[i]][data_t0[[i]]$gene_id == y, ][["stddev_norm"]], 
                                    Stddev_Norm)
        )
    }
  }
  # Enlist
  output <- list(data_s_s, data_s)
  names(output) <- c("data_s_s", "data_s")
  # Return
  return(output)
}


# ---------------------------------------------------------------------------- #
estimate_hl <- function(project, 
                        experiment, 
                        data_s, 
                        grps = NULL, 
                        time_limit = 16, 
                        warnings = FALSE) {
  # subset specifyed groups
  if (!is.null(grps)) {
    data_s <- data_s %>% filter(Group %in% grps)
  }
  # Set nls params
  nls_sets <- nls.control(maxiter = 300, 
                          tol = 1e-05, 
                          minFactor = 1/2048,
                          printEval = FALSE, 
                          warnOnly = warnings)
  # Set inner function
  fit_expdec <- function(grp, data_s) {
    ## define 
    iter <- 0
    trashcan <- tibble()
    data_p <- tibble()
    data_p_t <- tibble()
    ## subset group
    data_s_tmp <- data_s[data_s$Group == grp, ]
    # set up progress bar
    g_set <- length(unique(data_s_tmp$gene_id))
    pb <- txtProgressBar(min = 0, max = g_set, style = 3)
    iter <- 0
    for (y in unique(data_s_tmp$gene_id)) {
      iter <- iter + 1
      ### Estimate c_0 coeficient from data
      c_0 <- min(data_s_tmp[data_s_tmp$gene_id == y, ]$Conv_Rate_Norm) * 0.5 # based on https://stats.stackexchange.com/questions/160552/
      
      ### Estimate nls model coeficients from linearized model
      model_0 <- lm(log(Conv_Rate_Norm - c_0) ~ Timepoint, 
                    data = data_s_tmp[data_s_tmp$gene_id == y, ], 
                    na.action = "na.exclude")
      start <- list(a = exp(coef(model_0)[1]), 
                    b = coef(model_0)[2], 
                    c = c_0)
      
      ### Fit non-linear model
      model <- tryCatch(
        {
          nls(Conv_Rate_Norm ~ a * exp(b * Timepoint) + c, 
              data = data_s_tmp[data_s_tmp$gene_id == y, ], 
              start = start, 
              na.action = "na.exclude", 
              control = nls_sets)
        },
        error = function(e) {
          #message(paste0("\nExponential decay cannot be fitted for ", y, " in ", grp))
          #message("The original error message:")
          #message(e)
          return(NA)
        }
      )
      
      # skip
      if (length(model) == 1) {
        # collect trash
        trash <- tibble(gene_id = y)
        trashcan <- rbind(trashcan, trash)
        next
      }
      
      ### Predict values
      p.int <-  predict(model, 
                        data.frame(Timepoint = seq(0, time_limit, by = 0.01)),
                        interval = 'confidence',
                        level =0.95)
      
     
      
      
      dat_p <- tibble(gene_id = y,
                      gene_name = unique(data_s_tmp[data_s_tmp$gene_id == y, ]$gene_name), 
                      Group = grp,
                      Timepoint = seq(0, time_limit, by = 0.01), 
                      Pred_Conv_Rate_Norm = p.int[, 1],
                      Pred_Conv_Rate_Norm_UpInt = p.int[, 2],
                      Pred_Conv_Rate_Norm_LoInt = p.int[, 3])
      
      ### Estimate half-life
      dat_f <- dat_p %>% 
        mutate(Conv_Diff = (0.5 - Pred_Conv_Rate_Norm)^2) %>% 
        dplyr::filter(Conv_Diff == min(Conv_Diff)) %>% 
        dplyr::filter(Timepoint == min(Timepoint)) %>% # in case of extremely fast degradation
        dplyr::select(gene_id, gene_name, Group, Half_life = Timepoint) 
      ## & Error rates
      
      ### Bind
      data_p <- rbind(data_p, dat_f)
      data_p_t <- rbind(data_p_t, dat_p)
      ### Report
      #setTxtProgressBar(pb, iter)
      if (!as.logical(iter %% 250)) {
        message(iter, " transcripts processed in ", grp, "...\n")
      }
      if (iter == g_set) {
        message("Total ", iter, " transcripts processed")
      }
    }
    
    # Write aux results
    dummy <- tryCatch(
      {
        write_tsv(data_p_t, paste0("/local/artem/Projects/", project, "/Results/hl_estimation/", experiment, "/", grp, "_data_p_t.tsv"),
                  col_names = TRUE)
        write_tsv(trashcan, paste0("/local/artem/Projects/", project, "/Results/hl_estimation/", experiment, "/", grp, "_trashcan.tsv"),
                  col_names = TRUE)
      },
      error = function(cond) {
        message("Why though?")
        return(NA)
      }
    )
    rm(data_p_t)
    
    # Return main results
    return(data_p)
  }
  
  ## parallel over groups
  message("estimating half-lives...")
  tictoc::tic()
  library(doParallel)
  registerDoParallel(cores = length(unique(data_s$Group)))
  foreach(i = unique(data_s$Group)) %dopar% {
    message("starting ", i, ".....")
    fit_expdec(grp = i, data_s)
  } -> hl_l
  hl_t <- purrr::reduce(hl_l, rbind)
  tictoc::toc()
  
  
  # filter junk
  # for (z in 1:length(hl_l)) {
  #   message(paste0(length(unique(hl_l[[z]][["trash"]]$gene_id)), " genes have been removed from ", names(hl_l[z])))
  #   hl_l[[z]] <- hl_l[[z]][["hl"]] %>% dplyr::filter(!gene_id %in% hl_l[[z]][["trash"]]$gene_id)
  # }
  # hl <- purrr::reduce(hl_l, rbind)
  
  
  # Output
  return(hl_t)
}


