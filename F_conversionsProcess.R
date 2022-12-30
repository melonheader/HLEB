#!/usr/bin/env Rscript
### INFO: Auxiliary functions for conversionsProcess.R
### DATE: 20.08.2022
### AUTHOR: Artem Baranovskii

library(tidyverse)
library(magrittr)
library(janitor)


# Read in conversions at T0 for quality control
# ------------------------------------------------------------------------------
expand_grid_unique <- function(x, y, include.diagonal=FALSE) {
  x <- unique(x)
  y <- unique(y)
  g <- function(i) {
    z <- setdiff(y, x[seq_len(i - include.diagonal)])
    if(length(z)) cbind(x[i], z, deparse.level = 0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

# Read in conversions at T0 for quality control
# ------------------------------------------------------------------------------
comp_T0NN <- function(tab_comb) {
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
                    t.base <- data.table::fread(y, data.table = F) #%>% filter(V2 > base_cutoff) 
                    colnames(t.base) <- c("gene", "base_count")
                    t.conv <- data.table::fread(z, data.table = F) #%>% filter(V2 > conv_cutoff)
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
  ## Summarise genes
  t_out %>% 
    group_by(sample_name, group, timepoint, b.rep, conv, gene, base_count) %>% 
    summarise(conv_count = sum(conv_count), 
              .groups = "drop") -> t_out
  return(t_out)
}

# Plot conversions at T0 for quality control
# ------------------------------------------------------------------------------
plot_TONN.qc <- function(tab_cr_T0NN, percent_limit = 6, save_pdf = F, path_pdf = "") {
  
  ## keep genes that are present in all replicates
  gns.int <- purrr::reduce(split(tab_cr_T0NN$gene, tab_cr_T0NN$sample_name), intersect)
  tab_cr_T0NN <- tab_cr_T0NN %>% 
    filter(gene %in% gns.int) %>% 
    mutate(col_var = ifelse(conv == "TC", "red", "black")) #%>% 
    #filter((conv_count / base_count) * 100)
  ## plot
  ggplot(data = tab_cr_T0NN, 
         aes(x = conv, y = (conv_count / base_count) * 100, color = col_var)) + 
    geom_boxplot(outlier.shape = NA) + 
    scale_color_manual(values = c("grey30", "red"), guide = "none") +
    stat_summary(fun.data = function(x) tibble(y = -(0.25 * percent_limit), 
                                               label = round(mean(x), 2)), 
                 geom = "text", 
                 vjust = 0.5, 
                 hjust = 0.1, 
                 angle = 90) +
    #geom_hline(yintercept = 0) +
    coord_cartesian(ylim = c(-2, percent_limit)) +
    facet_grid(group ~ b.rep) +
    #ggforce::facet_zoom(ylim = c(0, 0.6), zoom.size = , zoom.data = ifelse(conversion != "T->C", NA, FALSE), horizontal = FALSE) +
    theme_bw(base_size = 14) + 
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust = 1, 
                                     colour = c(rep("black", 10), "red", "black"))
    ) +
    labs(y = "Conversion rate (%)", 
         x = "", subtitle = paste0("n genes = ", length(gns.int))) -> p.out
  
  if (save_pdf) {
    if (path_pdf == "") {
      stop("specify path in 'path_pdf' argument")
    }
    ggsave(plot = p.out,
           filename = file.path(path_pdf, "CR_T0NN.pdf"), 
           device = "pdf", width = 6, height = 6)
  } else {
    return(p.out)
  }
}

# Read in T -> C conversions 
# ------------------------------------------------------------------------------
comp_TC <- function(tab_comb) {
  tab_comb <- tab_comb %>% filter(str_detect(conv, "TC"))
  purrr::pmap_dfr(list(tab_comb$stem,
                       tab_comb$path_base,
                       tab_comb$path_conv,
                       tab_comb$path_snp,
                       tab_comb$conv), 
                  function(x, y, z, w, r) {
                    sample_name = str_split_fixed(x, "_trimmed", 2) %>% `[`(, 1)
                    grp = str_split_fixed(sample_name, "_", 2) %>% `[`(, 1)
                    tpoint = str_split_fixed(sample_name, "_", 3) %>% `[`(, 2)
                    br = str_split_fixed(sample_name, "_", 3) %>% `[`(, 3)
                    t.base <- data.table::fread(y, data.table = F) #%>% filter(V2 > base_cutoff)
                    colnames(t.base) <- c("gene", "base_count")
                    t.conv <- data.table::fread(z, data.table = F) #%>% filter(V2 > conv_cutoff)
                    colnames(t.conv) <- c("gene", "conv_count")
                    t.snp <- data.table::fread(w, data.table = F)
                    colnames(t.snp) <- c("gene", "snp_count")
                    t.merged <- inner_join(t.base, t.conv, by = "gene") %>% 
                      left_join(., t.snp, by = "gene") %>% 
                      mutate(conv = r,
                             sample_name = sample_name, 
                             group = grp, 
                             timepoint = tpoint,
                             b.rep = br,
                             snp_count = ifelse(is.na(snp_count), 0, snp_count)) %>% 
                      dplyr::select(sample_name, group, timepoint, b.rep, conv, everything())
                  }
  ) -> t_out
  ## Summarise genes
  t_out %>% 
    group_by(sample_name, group, timepoint, b.rep, conv, gene, base_count) %>% 
    summarise(conv_count = sum(conv_count), 
              snp_count = unique(snp_count), 
              .groups = "drop") -> t_out
  
  return(t_out)
}


# Filter T -> C conversions 
# ------------------------------------------------------------------------------
filt_TC <- function(tab_cr_TC, base_cutoff = 14, conv_cutoff = 5, 
                    min_T0 = 3, 
                    min_T4 = 2, 
                    min_T8 = 0,
                    min_T16 = 0, 
                    group_filt_action = "intersect",
                    return_stats = T, out_drop = NA) {
  ## record dropped entries
  trashcan <- list()
  trashcan$total <- dim(tab_cr_TC)[1]
  
  ## Remvoe outliers of specified
  if (!is.na(out_drop[1])) {
    tab_cr_TC <- tab_cr_TC %>% filter(!sample_name %in% out_drop)
  }
  
  ## Drop entries with more SNP counts than T -> C
  trashcan$snp <- sum(tab_cr_TC$snp_count > tab_cr_TC$snp_count)
  tab_cr_TC <- tab_cr_TC %>% filter(conv_count >= snp_count)
  ## Drop entries based on T -> C and T coverage thresholds
  trashcan$TC <- sum((tab_cr_TC$conv_count < conv_cutoff) | (tab_cr_TC$base_count < base_cutoff)) 
  tab_cr_TC <- tab_cr_TC %>% filter(conv_count >= conv_cutoff, 
                                    base_count >= base_cutoff)
  ## Drop genes based on detection at timepoints
  tab_cr_TC %>% 
    group_by(group, timepoint, conv, gene) %>% 
    summarize(n_obs = n(), .groups = "drop") %>% 
    pivot_wider(., names_from = "timepoint", values_from = "n_obs", values_fill = 0) %>% 
    dplyr::select(group, conv, gene, T0, T4, T8, T16) -> t.obs_TC
  gns.tokeep_Ts <- t.obs_TC %>% 
    filter(T0 >= min_T0, 
           T4 >= min_T4,
           T8 >= min_T8, 
           T16 >= min_T16) %$% 
    split(.$gene, .$group)
  if (group_filt_action == "intersect") {
    gns.int_Ts <- purrr::reduce(gns.tokeep_Ts, intersect)
    trashcan$obs_T <- t.obs_TC %>% mutate(dropped = !(gene %in% gns.int_Ts))
    tab_cr_TC <- tab_cr_TC %>% 
      filter(gene %in% gns.int_Ts) %>% 
      mutate(cr_perc = ((conv_count - snp_count) / base_count) * 100)
  } else {
    purrr::map2_dfr(names(gns.tokeep_Ts), 
                    gns.tokeep_Ts, 
                    ~ tab_cr_TC %>% 
                      filter(group == .x, gene %in% .y) %>% 
                      mutate(cr_perc = ((conv_count - snp_count) / base_count) * 100)
                    ) -> tab_cr_TC
    purrr::map2_dfr(names(gns.tokeep_Ts), 
                    gns.tokeep_Ts, 
                    ~ t.obs_TC %>% 
                      filter(group == .x) %>% 
                      mutate(dropped = !(gene %in% .y))
    ) -> trashcan$obs_T
  }
  if (return_stats) {
    list(
      tab_cr_TC.filt = tab_cr_TC,
      filt_stats = trashcan
    ) -> l.out
    return(l.out)
  } else {
    return(t.obs_TC)
  }
}


# Normalise T -> C conversions to timepoint zero
# ------------------------------------------------------------------------------
norm_TC <- function(tab_cr_TC_filtered) {
  ## estimate T0 averages
  tab_cr_TC_filtered %>% 
    filter(timepoint == "T0") %>% 
    group_by(group, gene) %>% 
    summarise(cr_perc.Avg_T0 = mean(cr_perc), 
              .groups = "drop") -> T0_avgs
  tab_cr_TC_filtered %>% 
    filter(timepoint != "T0") %>% 
    dplyr::select(group, timepoint, b.rep, gene, cr_perc) %>% 
    pivot_wider(., names_from = c("timepoint", "b.rep"),
                values_from = cr_perc, names_glue = "{timepoint}_{b.rep}") %>% 
    dplyr::select(group, gene, contains("T0"), contains("T4"), contains("T8"), contains("T16")) -> t.cr_wide
  ##
  apply(t.cr_wide %>% 
          dplyr::select(-group, -gene),
        2,
        function(x) x / T0_avgs$cr_perc.Avg_T0) -> t.cr_wide.norm
  
  cbind(t.cr_wide[, 1:2], 
        t.cr_wide.norm) %>% 
    pivot_longer(., contains("T"), names_to = "tmp", values_to = "cr_norm") %>% 
    mutate(timepoint = str_split_fixed(tmp, "_", 2)[, 1],
           b.rep = str_split_fixed(tmp, "_", 2)[, 2]) %>% 
    dplyr::select(-tmp) -> t.cr_long.norm
  ##
  t.out <- tab_cr_TC_filtered %>%  
    dplyr::select(group, gene, timepoint, b.rep, cr_raw = cr_perc) %>% 
    left_join(., t.cr_long.norm, by = c("group", "gene", "timepoint", "b.rep")) %>% 
    mutate(cr_norm = ifelse(timepoint == "T0", 1, cr_norm),
           timepoint = as.numeric(str_remove(timepoint, "T")))
  return(t.out)
}


# Function to estimate mRNA Half-lives from T -> C conversion rates data
# ------------------------------------------------------------------------------
estimate_HL <- function(tab_cr_TC_norm, time_limit = 24, n_cores = 8, out_path) {
  ## define 
  require(purrr)
  require(furrr)
  ##
  if (!dir.exists(file.path(out_path, "hl_estimation"))) {
    dir.create(file.path(out_path, "hl_estimation"))
  } else {
    sapply(list.files(file.path(out_path, "hl_estimation"), full.names = T), file.remove)
  }
  ## report input data
  t.tr.stats <- tab_cr_TC_norm %>% 
    group_by(group) %>% 
    distinct(gene) %>% 
    janitor::tabyl(group) %>% as_tibble()
  message(paste0("Input data:\n"))
  message(paste0("Group: ", t.tr.stats$group, ";\t", "n unique transcripts = ", t.tr.stats$n, "\n"))
  
  
  ## create a dummy var for splitting
  tab_cr_TC_norm <- tab_cr_TC_norm %>% mutate(dummy = paste0(group, "=", gene))
  ## ## remove obs with cr_norm > 1 after T0
  gns.drop <- tab_cr_TC_norm %>% 
    filter(timepoint > 0, cr_norm > 1) %>% 
    distinct(dummy) %>% pull(dummy)
  tab_cr_TC_norm <- tab_cr_TC_norm %>% filter(!dummy %in% gns.drop)
  uniq.dummy <- unique(tab_cr_TC_norm$dummy)
  ## split data based on n cores provided
  iter_size = length(uniq.dummy) %/% n_cores
  iter_rmdr = length(uniq.dummy) %% n_cores
  v.splits =  c(rep(1:n_cores, each = iter_size), rep(n_cores, iter_rmdr))
  ##
  tab_cr_TC_norm <- tab_cr_TC_norm %>% inner_join(., tibble(dummy = uniq.dummy, dummy_1 = v.splits), by = "dummy")
  
  ## split data to parallelize
  l_cr_TC_norm <- split(tab_cr_TC_norm, tab_cr_TC_norm$dummy_1)
  
  ## set up workers
  future::plan(multisession, workers = n_cores)
  furrr::furrr_options(stdout = T, seed = T)
  furrr::future_map_dfr(l_cr_TC_norm, 
                        function(x) { fit_expdec(x, time_limit, out_path) }, .options = furrr_options(seed = T) 
                        ) -> t.out
  return(t.out)
}

# ------------------------------------------------------------------------------
fit_expdec <- function(split_cr_TC_norm, time_limit, out_path) {
  # generate random run code
  runcode <- paste0(sample(letters, 2), sample(1:9, 2), collapse = "")
  
  #
  if (!dir.exists(file.path(out_path, "hl_estimation"))) {
    dir.create(dir.exists(file.path(out_path, "hl_estimation")))
  }
  
  # Set nls params
  nls_sets <- nls.control(maxiter = 300, 
                          tol = 1e-05, 
                          minFactor = 1/2048,
                          printEval = FALSE, 
                          warnOnly = F)
  ## define 
  trashcan <- tibble()
  data_p <- tibble()
  data_p_t <- tibble()
  # set up progress bar
  uniq.obs <- unique(split_cr_TC_norm$dummy)
  pb <- txtProgressBar(min = 0, max = length(uniq.obs), style = 3)
  iter <- 0
  for (y in uniq.obs) {
    iter <- iter + 1
    ##
    grp = str_split_fixed(y, "=", 2)[, 1]
    gene = str_split_fixed(y, "=", 2)[, 2]
    split_tmp <- split_cr_TC_norm %>% 
      dplyr::filter(dummy == y) %>%
      dplyr::select(-contains("dummy"))
    
    ### Estimate c_0 coeficient from data
    c_0 <- min(split_tmp$cr_norm) * 0.5 # based on https://stats.stackexchange.com/questions/160552/
    
    ### Estimate nls model coeficients from linearized model
    model_0 <- lm(log(cr_norm - c_0) ~ timepoint, 
                  data = split_tmp, 
                  na.action = "na.exclude")
    start <- list(a = exp(coef(model_0)[1]), 
                  b = coef(model_0)[2], 
                  c = c_0)
    
    ### Fit non-linear model
    model <- tryCatch(
      {
        nls(cr_norm ~ a * exp(b * timepoint) + c, 
            data = split_tmp, 
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
      trash <- tibble(dummy = y, dropped_by = 1)
      trashcan <- rbind(trashcan, trash)
      next
    }
    ### Predict values
    p.int <- tryCatch(
      {
        investr::predFit(model, 
                         newdata = data.frame(timepoint = seq(0, time_limit, by = 0.025)), 
                         interval = "confidence",
                         level = 0.95)
      },
      error = function(e) {
        return(NA)
      }
    )
    # skip
    if (length(p.int) == 1) {
      # collect trash
      trash <- tibble(dummy = y, dropped_by = 2)
      trashcan <- rbind(trashcan, trash)
      next
    }
    
    dat_p <- tibble(group = grp,
                    gene = gene,
                    timepoint = seq(0, time_limit, by = 0.025), 
                    cr_pred_norm = p.int[, 1],
                    cr_pred_norm_LoInt = p.int[, 2],
                    cr_pred_norm_UpInt = p.int[, 3])
    ## we cannot reliably estiamte error for half-lives that reach above 24 hours
    ## my approach was inherently wrong -- use linear interpolation
    # y = a*e^(b*T) + c -> ln(y - c) = b*T + ln(a) -> T = (ln(y - c) - ln(a)) / b
    ## new fit
    slope <- function(x, y) {
      return(cov(x, y) / var(x))
    }
    intercept <- function(x, y, slope){
      b <- mean(y) - (slope * mean(x))
      return(b)
    }
    tibble(x = dat_p$timepoint,
           y = log(dat_p$cr_pred_norm - coef(model)[["c"]]), 
    ) -> tmp_lin
    tmp_slp <- slope(tmp_lin$x, tmp_lin$y)
    tmp_int <- intercept(tmp_lin$x, tmp_lin$y, slope(tmp_lin$x, tmp_lin$y))
    
    ## Compute quadratic difference 
    ph_val = (0.5 - dat_p$cr_pred_norm)^2
    ## And estimate half-life and it's intervals
    hl_val = unique(dat_p$timepoint[which(ph_val == min(ph_val))])[1]
    if (hl_val < 24) {
      if (dat_p[dat_p$timepoint == hl_val, ]$cr_pred_norm_LoInt < coef(model)[["c"]]) {
        hl_UpInt = time_limit
      } else {
        hl_UpInt = (log(dat_p[dat_p$timepoint == hl_val, ]$cr_pred_norm_LoInt - coef(model)[["c"]]) - tmp_int) / tmp_slp
      }
      hl_LoInt = (log(dat_p[dat_p$timepoint == hl_val, ]$cr_pred_norm_UpInt - coef(model)[["c"]]) - tmp_int) / tmp_slp
    } else {
      hl_LoInt = NA
      hl_UpInt = NA
    }
    ## Assemble into a reporting table
    dat_f <- tibble(group = grp,
                    gene = gene,
                    Half_life = hl_val,
                    HL_LoInt95 = hl_LoInt,
                    HL_UpInt95 = hl_UpInt)
    
    ### Bind
    data_p <- rbind(data_p, dat_f)
    data_p_t <- rbind(data_p_t, dat_p)
    
    ### Report
    setTxtProgressBar(pb, iter)
  }
  
  # Write aux results
  dummy <- tryCatch(
    {
      write_tsv(data_p, file.path(out_path, "hl_estimation",  paste0(runcode, ".data_hl.tsv")),
                append = file.exists(file.path(out_path, "hl_estimation",  paste0(runcode, ".data_hl.tsv"))))
      write_tsv(data_p_t, file.path(out_path, "hl_estimation",  paste0(runcode, ".data_cr_pred.tsv")),
                append = file.exists(file.path(out_path, "hl_estimation",  paste0(runcode, ".data_cr_pred.tsv"))))
      write_tsv(trashcan, file.path(out_path, "hl_estimation",  paste0(runcode, ".trashcan.tsv")),
                append = file.exists(file.path(out_path, "hl_estimation",  paste0(runcode, ".trashcan.tsv"))))
    },
    error = function(cond) {
      message("Why though?")
      return(NA)
    }
  )
  
  # Return main results
  return(data_p)
}


# Plot T -> C conversions as boxes per replicate per timepoint
# ------------------------------------------------------------------------------
plot_TC.qc_boxes <- function(tab_cr_TC, y_lims = c(0, 20), save_pdf = F, path_pdf = "") {
  ggplot(data = tab_cr_TC %>% mutate(timepoint = factor(str_remove(timepoint, "T"), 
                                                        levels = c("0", "4", "8", "16"))), 
         aes(y = ((conv_count - snp_count) / base_count) * 100,
             x = timepoint)) + 
    geom_boxplot(aes(fill = group, color = b.rep), position = "dodge", outlier.shape = NA) + 
    scale_color_manual(values = rep("gray30", length(unique(tab_cr_TC$b.rep)))) +
    guides(color = "none") +
    coord_cartesian(ylim = y_lims) +
    labs(y = "T -> C conversion rate, %", x = "Timpoint, hours") +
    theme_bw(base_size = 14) -> p.out
  if (save_pdf) {
    if (path_pdf == "") {
      stop("specify path in 'path_pdf' argument")
    }
    ggsave(plot = p.out,
           filename = file.path(path_pdf, "CR_T0NN.pdf"), 
           device = "pdf", width = 6, height = 6)
  } else {
    return(p.out)
  }
}

# Plot KS test D statistic that compares T -> C CR between each pair of replicates
# ------------------------------------------------------------------------------
plot_TC.qc_ks <- function(tab_cr_TC, save_pdf = F, path_pdf = "") {
  
  tab_cr_TC <- tab_cr_TC %>% mutate(cr.snpCor_perc = ((conv_count - snp_count) / base_count) * 100)
  
  ## test distributions against each other within group, within timepoint
  br.names <- unique(tab_cr_TC$b.rep)
  br.grid <- expand_grid_unique(br.names, br.names, include.diagonal = T) %>% 
    as_tibble(., .name_repair = ~ make.names(., unique = T))
  colnames(br.grid) <- c("X", "Y")
  ##
  purrr::map_dfr(split(tab_cr_TC, tab_cr_TC$group),
                 function(split_gr) {
                   purrr::map_dfr(split(split_gr, split_gr$timepoint),
                                  function(split_gr_tp) {
                                    purrr::map2_dfr(br.grid$X,
                                                    br.grid$Y,
                                                    function(x, y) {
                                                      if (dim(split_gr_tp[split_gr_tp$b.rep == x, ])[1] == 0 |
                                                          dim(split_gr_tp[split_gr_tp$b.rep == y, ])[1] == 0) {
                                                        tibble(group = unique(split_gr_tp$group),
                                                               timepoint = unique(split_gr_tp$timepoint),
                                                               x = x,
                                                               y = y,
                                                               test = "KS",
                                                               stat = NA,
                                                               p.val = NA)
                                                      } else {
                                                        tmp.ks <- ks.test(split_gr_tp[split_gr_tp$b.rep == x, ]$cr.snpCor_perc,
                                                                          split_gr_tp[split_gr_tp$b.rep == y, ]$cr.snpCor_perc)
                                                        tibble(group = unique(split_gr_tp$group),
                                                               timepoint = unique(split_gr_tp$timepoint),
                                                               x = x,
                                                               y = y,
                                                               test = "KS",
                                                               stat = tmp.ks$statistic,
                                                               p.val = ifelse(tmp.ks$p.value == 0, 2.2e-16, tmp.ks$p.value))
                                                      }
                                                      
                                                    })
                                  })
                 }) %>% mutate(p.adj = p.adjust(p.val, method = "bonferroni"),
                               timepoint = factor(timepoint,
                                                  levels = c("T0", "T4", "T8", "T16")
                               )
                 ) -> tmp.p
  
  ggplot(data = tmp.p %>% mutate(timepoint = factor(str_remove(timepoint, "T"), 
                                                    levels = c("0", "4", "8", "16"))), 
         aes(x, y, fill = stat)) + 
    geom_tile(color = "black", 
              width = 0.96, height = 0.96) +
    scale_fill_gradient(name = "D statistic", low = "palegoldenrod", high = "orangered2", na.value = NA) + 
    facet_grid(group ~ timepoint) + 
    labs(x = "", y = "") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> p.out
  if (save_pdf) {
    if (path_pdf == "") {
      stop("specify path in 'path_pdf' argument")
    }
    ggsave(plot = p.out,
           filename = file.path(path_pdf, "CR_T0NN.pdf"), 
           device = "pdf", width = 6, height = 6)
  } else {
    return(p.out)
  }
}
