## Funcitons
# ---------------------------------------------------------------------------- #
estimate_hl <- function(project, 
                        experiment, 
                        data_s, 
                        grps = NULL, 
                        time_limit = 16, 
                        warnings = FALSE,
                        aux_results = FALSE) {
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
    data_p_t.filt <- tibble()
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
      p.int <- investr::predFit(model, 
                                newdata = data.frame(Timepoint = seq(0, time_limit, by = 0.025)), 
                                interval = "confidence",
                                level= 0.95)
      
      dat_p <- tibble(gene_id = y,
                      gene_name = unique(data_s_tmp[data_s_tmp$gene_id == y, ]$gene_name), 
                      Group = grp,
                      Timepoint = seq(0, time_limit, by = 0.025), 
                      Pred_Conv_Rate_Norm = p.int[, 1],
                      Pred_Conv_Rate_Norm_LoInt = p.int[, 2],
                      Pred_Conv_Rate_Norm_UpInt = p.int[, 3])
      
      ## we cannot reliably estiamte error for half-lives that reach above 24 hours
      ## use linear interpolation for error estimation
      # y = a*e^(b*T) + c -> ln(y - c) = b*T + ln(a) -> T = (ln(y - c) - ln(a)) / b
      ## new fit
      slope <- function(x, y) {
        return(cov(x, y) / var(x))
      }
      intercept <- function(x, y, slope){
        b <- mean(y) - (slope * mean(x))
        return(b)
      }
      tibble(x = dat_p$Timepoint,
             y = log(dat_p$Pred_Conv_Rate_Norm - coef(model)[["c"]]), 
             ) -> tmp_lin
      tmp_slp <- slope(tmp_lin$x, tmp_lin$y)
      tmp_int <- intercept(tmp_lin$x, tmp_lin$y, slope(tmp_lin$x, tmp_lin$y))
      
      
      ## Compute quadratic difference 
      ph_val = (0.5 - dat_p$Pred_Conv_Rate_Norm)^2
      ## And estimate half-life and it's intervals
      hl_val = dat_p$Timepoint[which(ph_val == min(ph_val))]
      if (hl_val < 24) {
        hl_LoInt = (log(dat_p[dat_p$Timepoint == hl_val, ]$Pred_Conv_Rate_Norm_UpInt - coef(model)[["c"]]) - tmp_int) / tmp_slp
        hl_UpInt = (log(dat_p[dat_p$Timepoint == hl_val, ]$Pred_Conv_Rate_Norm_LoInt - coef(model)[["c"]]) - tmp_int) / tmp_slp
      } else {
        hl_LoInt = NA
        hl_UpInt = NA
      }
      dat_p.filt <- dat_p %>% filter(Timepoint == hl_val)
      
      
      
      ## Assemble into a reporting table
      dat_f <- tibble(gene_id = y,
                      gene_name = unique(data_s_tmp[data_s_tmp$gene_id == y, ]$gene_name),
                      Group = grp,
                      Half_life = hl_val,
                      HL_LoInt95 = hl_LoInt,
                      HL_UpInt95 = hl_UpInt)
      
      ### Bind
      data_p <- rbind(data_p, dat_f)
      data_p_t <- rbind(data_p_t, dat_p)
      data_p_t.filt <- rbind(data_p_t.filt, dat_p.filt)
      
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
    if (aux_results) {
      dummy <- tryCatch(
        {
          write_tsv(data_p_t, paste0("/local/artem/Projects/", project, "/Results/hl_estimation/", experiment, "/", grp, "_data_p_t.tsv"),
                    col_names = TRUE)
          write_tsv(data_p_t, paste0("/local/artem/Projects/", project, "/Results/hl_estimation/", experiment, "/", grp, "_data_p_t.filt.tsv"),
                    col_names = TRUE)
          write_tsv(trashcan, paste0("/local/artem/Projects/", project, "/Results/hl_estimation/", experiment, "/", grp, "_trashcan.tsv"),
                    col_names = TRUE)
        },
        error = function(cond) {
          message("Why though?")
          return(NA)
        }
      )
    }
    rm(data_p_t)
    
    # Return main results
    return(data_p)
  }
  
  ## parallel over groups --- parallelize over transcripts, not groups
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
  
  
  filter junk
  for (z in 1:length(hl_l)) {
    message(paste0(length(unique(hl_l[[z]][["trash"]]$gene_id)), " genes have been removed from ", names(hl_l[z])))
    hl_l[[z]] <- hl_l[[z]][["hl"]] %>% dplyr::filter(!gene_id %in% hl_l[[z]][["trash"]]$gene_id)
  }
  hl <- purrr::reduce(hl_l, rbind)
  
  
  # Output
  return(hl_t)
}