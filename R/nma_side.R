#' Side-splitting model for network meta-analysis
#' 
#' @description
#' Side-splitting model for network meta-analysis
#' 
#' @param s Column name for study (default "s")
#' @param y Column name for reported outcome (default "y")
#' @param v Column name for outcome variance (default "v")
#' @param trt Column name for treatment (default "trt"); s, y, v, trt
#'   can be input as vectors directly if data is set to NULL
#' @param comp Length-2 vector listing the treatments in the contrast
#'   of interest; effect is calculated as the second treatment - the
#'   first treatment
#' @param ref Takes values "1", "2" or "sym"; "1" and "2" sets
#'   side-splitting models with comp[1] or comp[2] as reference
#'   treatment. "sym" sets the symmetric side-splitting model (default
#'   "sym")
#' @param data Input data frame (long form, i.e. one study-treatment
#'   pair per row)
#' @param printout Boolean for printing out the results (default TRUE)
#' @param large_var Large variance for incidental parameter problem
#'   (default 1e5)
#' @param conf confidence level of the confidence intervals (default
#'   .95); confidence intervals are constructed via likelihood-based
#'   approach
#' 
#' @details
#' Side-splitting model for network meta-analysis
#'
#' !! Currently only support continuous outcomes !!
#'
#' !! Check inconsistency identifiablity (van Valkenhoef, 2012) beforehand!!
#' 
#' @return
#' List ...
#' 
#' @author Ming-Chieh Shih \email{mcshih@mx.nthu.edu.tw}
#' 
#' @seealso \code{\link{nma}}, \code{\link{nma_incon}}
#' 
#' @references
#' Shih ...
#' Title.
#' \emph{Journal},
#' \bold{x}, p1--p2
#' 
#' @examples
#' library("dplyr")
#' library("tidyr")
#' 
#' #### Dong2013 LABA-ICS compared to LABA ####
#' 
#' data(Dong2013, package = "netmeta")
#' 
#' Dong_new <- as_tibble(Dong2013) %>% 
#'   group_by(id) %>% 
#'   summarize(have_zero = (prod(death) == 0) | (prod(randomized - death) ==0), 
#'             all_zero = (sum(death) == 0) | (sum(randomized - death) ==0)) %>% 
#'   filter(!all_zero) %>%
#'   inner_join(Dong2013) %>%
#'   mutate(y = log((death + 0.5 * have_zero) / (randomized - death + 0.5 * have_zero)),
#'          v = 1/(death + 0.5 * have_zero) + 1/(randomized - death + 0.5 * have_zero))
#' 
#' Dong_res_sym <- nma_side(y = "y", v = "v", trt = "treatment", s = "id", 
#'                          comp = c("LABA", "LABA-ICS"), 
#'                          data = Dong_new)
#' Dong_res_laba <- nma_side(y = "y", v = "v", trt = "treatment", s = "id", 
#'                           comp = c("LABA", "LABA-ICS"), ref = "1", 
#'                           data = Dong_new)
#' Dong_res_laba_ics <- nma_side(y = "y", v = "v", trt = "treatment", s = "id", 
#'                               comp = c("LABA", "LABA-ICS"), ref = "2", 
#'                               data = Dong_new)
#' 
#' #### Linde2015 Hypericum compared to Placebo ####
#' 
#' data(Linde2015, package = "netmeta")
#' 
#' Linde_new <- as_tibble(Linde2015) %>%
#'   select(id, starts_with("treatment"), starts_with("resp"), starts_with("n")) %>%
#'   pivot_longer(treatment1:n3, names_to = ".value", names_pattern = "(.*).") %>%
#'   filter(!is.na(resp))
#' 
#' Linde_new <- Linde_new %>% 
#'   group_by(id) %>% 
#'   summarize(have_zero = (prod(resp) == 0) | (prod(resp - n) ==0), 
#'             all_zero = (sum(resp) == 0) | (sum(resp - n) ==0)) %>% 
#'   filter(!all_zero) %>%
#'   inner_join(Linde_new) %>%
#'   mutate(y = log((resp + 0.5 * have_zero) / (n - resp + 0.5 * have_zero)),
#'          v = 1/(resp + 0.5 * have_zero) + 1/(n - resp + 0.5 * have_zero))
#' 
#' Linde_res_sym <- nma_side(y = "y", v = "v", trt = "treatment", s = "id", 
#'                           comp = c("Placebo", "Hypericum"), 
#'                           data = Linde_new)
#' Linde_res_placebo <- nma_side(y = "y", v = "v", trt = "treatment", s = "id", 
#'                               comp = c("Placebo", "Hypericum"), ref = "1", 
#'                               data = Linde_new)
#' Linde_res_hypericum <- nma_side(y = "y", v = "v", trt = "treatment", s = "id", 
#'                                 comp = c("Placebo", "Hypericum"), ref = "2", 
#'                                 data = Linde_new)
#' 
#' ### Pagliaro1992 (the data used in Shih & Tu, 2021) ###
#' ### Beta-blocker compared to Sclerotherapy          ###
#' 
#' data(Pagliaro1992 , package = "meta")
#' Pagliaro_new <- as_tibble(Pagliaro1992) %>% 
#'   group_by(id) %>% 
#'   slice(1) %>% 
#'   ungroup() %>%
#'   select(id, bleed.plac, n.plac) %>% 
#'   rename(bleed.exp = bleed.plac, n.exp = n.plac) %>%
#'   mutate(treat.exp = "Placebo") %>%
#'   bind_rows(as_tibble(Pagliaro1992) %>% 
#'               select(id, treat.exp, bleed.exp, n.exp))
#' 
#' Pagliaro_new <- Pagliaro_new %>% 
#'   group_by(id) %>% 
#'   summarize(have_zero = (prod(bleed.exp) == 0) | (prod(n.exp - bleed.exp) ==0), 
#'             all_zero = (sum(bleed.exp) == 0) | (sum(n.exp - bleed.exp) ==0)) %>% 
#'   filter(!all_zero) %>%
#'   inner_join(Pagliaro_new) %>%
#'   mutate(y = log((bleed.exp + 0.5 * have_zero) / (n.exp - bleed.exp + 0.5 * have_zero)),
#'          v = 1/(bleed.exp + 0.5 * have_zero) + 1/(n.exp - bleed.exp + 0.5 * have_zero))
#' 
#' Pagliaro_res_sym <- nma_side(y = "y", v = "v", trt = "treat.exp", s = "id", 
#'                              comp = c("Sclerotherapy", "Beta-blocker"), 
#'                              data = Pagliaro_new)
#' Pagliaro_res_sclero <- nma_side(y = "y", v = "v", trt = "treat.exp", s = "id", 
#'                                 comp = c("Sclerotherapy", "Beta-blocker"), ref = "1",
#'                                 data = Pagliaro_new)
#' Pagliaro_res_beta <- nma_side(y = "y", v = "v", trt = "treat.exp", s = "id", 
#'                               comp = c("Sclerotherapy", "Beta-blocker"), ref = "2",
#'                               data = Pagliaro_new)
#' 
#' @export nma_side


nma_side <- function(s = "s", y = "y", v = "v", trt = "trt", comp, ref = "sym", 
                      data = NULL, printout = TRUE, large_var = 1e5, conf = .95){
  if(is.null(data)) data <- sys.frame(sys.parent())
  if(is.numeric(ref)) ref <- as.character(ref)
  
  ## Get rid of warnings 'no visible binding for global variable'
  . <- dmat <- incmat <- sigmat <- NULL
  
  # Obtain coding book for treatment recoding indices and number of treatments
  trt_list <- levels(factor(data[[trt]]))
  ntrt <- length(trt_list)
  
  # Obtain treatment indices for the contrast of interest (coi)
  coi <- c(which(trt_list == comp[1]), which(trt_list == comp[2]))
  
  # coi contrast vector for later use
  cvec <- rep(0, ntrt); cvec[coi[1]] <- -1; cvec[coi[2]] <- 1
  
  # Pivot to wide table, define direct study with coi, and 
  # fill in placeholder variance for missing arms
  data <- data.frame(s = data[[s]], 
                     y = data[[y]], 
                     v = data[[v]], 
                     trt = as.numeric(factor(data[[trt]])))
  data <- data %>% 
    pivot_wider(id_cols = s, names_from = trt, 
                names_sep = "", values_from = c(y, v)) %>%
    mutate(., dir = as.numeric((!is.na(.[,paste0("y", coi[1])])) & 
                                 (!is.na(.[,paste0("y", coi[2])])))) %>%
    mutate(across(starts_with('v'), ~replace_na(., 1)))
  
  # Variable name vectors for future use
  y_var <- paste0("y", 1:ntrt)
  u_var <- paste0("u", 1:ntrt)
  v_var <- paste0("data.v", 1:ntrt)
  
  # Define dat for direct and indirect studies
  data_dir <- mxData(data.frame(data %>% filter(dir == 1)), type = "raw")
  data_ind <- mxData(data.frame(data %>% filter(dir == 0)), type = "raw")
  
  # Definition of parameter matrices
  # d (d1 - d'ntrt'): column vector of size ntrt for treatment effects 
  # incmat (inc): inconsistency parameter
  # sigsq (sigsq): half the heterogeneity variance
  d <- mxMatrix(type = "Full", nrow = ntrt, ncol = 1, free = c(F, rep(T, ntrt-1)),
                labels = c(paste0("d", 1:ntrt)), values = rep(0, ntrt), name = "dmat")
  inc <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, 
                  labels = "inc", name = "incmat", values = 0)
  sig <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, labels="sig", name = "sigmat")
  sigsq <- mxAlgebra(sigmat^2, name = "sigsq")
  
  # Model lines for study effects
  s_v <- mxPath(from = "mu", arrows = 2, free = FALSE, values = large_var)
  s_p <- mxPath(from = "mu", to = y_var, free = FALSE, values = 1)
  
  # Model lines for treatment effects
  t_p <- mxPath(from = "one", to = y_var, free = c(F, rep(T, ntrt-1)), 
                labels = paste0("d", 1:ntrt))
  
  # Model lines for heterogeneity effects
  h_v <- mxPath(from = u_var, arrows = 2, free = FALSE, labels = "sigsq[1,1]")
  h_p <- mxPath(from = u_var, to = y_var, free = FALSE, values = 1)
  
  # Model lines for error variance
  e_v <- mxPath(from = y_var, arrows = 2, free = FALSE, labels = v_var)
  
  # Model lines for inconsistency 
  w_v <- mxPath(from = "w", arrows = 2, free = FALSE, values = 0)
  w_m <- mxPath(from = "one", to = "w", free = TRUE, labels = "inc")
  
  if(ref == "1"){
    w_p <- mxPath(from = "w", to = paste0("y", coi), free = FALSE, values = c(0, 1))
  }else if(ref == "2"){
    w_p <- mxPath(from = "w", to = paste0("y", coi), free = FALSE, values = c(-1, 0))
  }else if(ref == "sym"){
    w_p <- mxPath(from = "w", to = paste0("y", coi), free = FALSE, values = c(-0.5, 0.5))
  }
  
  # Additional calculations for direct and indirect estimates
  b <- mxMatrix(type = "Full", nrow = 1, ncol = ntrt, free = FALSE,
                values = cvec, name = "b")
  dir_est <- mxAlgebra(b %*% dmat + incmat, name = "dir_est")
  ind_est <- mxAlgebra(b %*% dmat, name = "ind_est")
  
  # Putting together the models for direct and indirect studies
  mod_dir <- mxModel("dir", type="RAM", 
                     manifestVars = y_var, latentVars = c(u_var, "mu", "w"),
                     data_dir,
                     s_v, s_p, t_p, h_v, h_p, e_v, d, sig, sigsq, 
                     inc, w_v, w_m, w_p)
  mod_ind <- mxModel("ind", type="RAM", 
                     manifestVars = y_var, latentVars = c(u_var, "mu"),
                     data_ind,
                     s_v, s_p, t_p, h_v, h_p, e_v, d, sig, sigsq)
  mod <- mxModel("mod", mod_dir, mod_ind, 
                 mxFitFunctionMultigroup(c("dir.fitfunction","ind.fitfunction")),
                 b, d, inc, dir_est, ind_est, 
                 mxCI(c('dir_est','ind_est','inc'), interval = conf))
  
  res <- mxRun(mod, intervals = TRUE)
  output <- res$output$confidenceIntervals[, c(2,1,3)]
  rownames(output) <- c("Direct Estimate", "Indirect Estimate", "Difference")
  colnames(output) <- c("Estimate", "Lower Bound", "Upper Bound")
  if(printout){
    if(ref == "1"){
      cat("Symmetric side-splitting model for ", comp[2], " compared to ", 
          comp[1], ", ", comp[1], " as reference:\n", sep = "")
    }else if(ref == "2"){
      cat("Symmetric side-splitting model for ", comp[2], " compared to ", 
          comp[1], ", ", comp[1], " as reference:\n", sep = "")
    }else if(ref == "sym"){
      cat("Symmetric side-splitting model for ", comp[2], " compared to ",
          comp[1], ":\n", sep = "")
    }
    print(output)
  }
  return(res)
}
