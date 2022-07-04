# =========================================================================== #
#               Evidence-splitting model for network meta-analysis            #
# --------------------------------------------------------------------------- #
#             !! Currently only support continous outcomes !!                 #
#  !! Check inconsistency identifiablity (van Valkenhoef, 2012) beforehand!!  #
#                                                                             #
#  main function: nma_incon(s, y, v, trt, comp, data, printout, large_var)    #
#                                                                             #
#  data: Input data frame (long form, i.e. one study-treatment pair per row)  #
#  s: Column name for study (default "s")                                     #
#  y: Column name for reported outcome (default "y")                          #
#  v: Column name for outcome variance (default "v")                          #
#  trt: Column name for treatment (default "trt")                             #
#    -- s, y, v, trt can be input as vectors directly if data is set to NULL  #
#  comp: Length-2 vector listing the treatments in the contrast of interest   #
#    -- Effect is calculated as the second treatment - the first treatment    #
#  printout: Boolean for printing out the results (default T)                 #
#  large_var: Large variance for incidental parameter problem (default 1e5)   #
#  conf: confidence level of the confidence intervals (default .95)           #
#    -- Confidence intervals are constructed via likelihood-based approach    #
#                                                                             #
# --------------------------------------------------------------------------- #
#  Author : Ming-Chieh Shih (Last updated 2022-07-03)                         #
# =========================================================================== #

require.i <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# netmeta and meta are called for real-world data 
# (Dong2013, Linde 2015, Pagliaro1992)
require.i(c("OpenMx", "tidyverse", "netmeta", "meta"))

nma_incon <- function(s = "s", y = "y", v = "v", trt = "trt", comp, 
                      data = NULL, printout = T, large_var = 1e5, conf = .95){
  if(is.null(data)) data <- sys.frame(sys.parent())
  
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
  inc <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = T, 
                  labels = "inc", name = "incmat", values = 0)
  sig <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = T, labels="sig", name = "sigmat")
  sigsq <- mxAlgebra(sigmat^2, name = "sigsq")
  
  # Model lines for study effects
  s_v <- mxPath(from = "mu", arrows = 2, free = F, values = large_var)
  s_p <- mxPath(from = "mu", to = y_var, free = F, values = 1)
  
  # Model lines for treatment effects
  t_p <- mxPath(from = "one", to = y_var, free = c(F, rep(T, ntrt-1)), 
                labels = paste0("d", 1:ntrt))
  
  # Model lines for heterogeneity effects
  h_v <- mxPath(from = u_var, arrows = 2, free = F, labels = "sigsq[1,1]")
  h_p <- mxPath(from = u_var, to = y_var, free = F, values = 1)
  
  # Model lines for error variance
  e_v <- mxPath(from = y_var, arrows = 2, free = F, labels = v_var)
  
  # p: vector of inconsistency parameter splitting weight 
  #    between the treatments of interest
  v <- mxMatrix(type = "Full", nrow = 2, ncol = 1, free = F, 
                labels = paste0("data.v", coi), name = "v")
  a <- mxMatrix(type = "Full", nrow = 2, ncol = 1, free = F, 
                values = c(-1, 1), name = "a")
  p <- mxAlgebra(a * (v + sigsq) / (sum(v) + 2*sigsq), name = "p")
  
  # Model lines for inconsistency 
  w_v <- mxPath(from = "w", arrows = 2, free = F, values = 0)
  w_m <- mxPath(from = "one", to = "w", free = T, labels = "inc")
  w_p <- mxPath(from = "w", to = paste0("y", coi), free = F, labels = c("p[1,1]", "p[2,1]"))
  
  # Additional calculations for direct and indirect estimates
  b <- mxMatrix(type = "Full", nrow = 1, ncol = ntrt, free = F,
                values = cvec, name = "b")
  dir_est <- mxAlgebra(b %*% dmat + incmat, name = "dir_est")
  ind_est <- mxAlgebra(b %*% dmat, name = "ind_est")
  
  # Putting together the models for direct and indirect studies
  mod_dir <- mxModel("dir", type="RAM", 
                     manifestVars = y_var, latentVars = c(u_var, "mu", "w"),
                     data_dir,
                     s_v, s_p, t_p, h_v, h_p, e_v, d, sig, sigsq, 
                     inc, v, a, p, w_v, w_m, w_p)
  mod_ind <- mxModel("ind", type="RAM", 
                     manifestVars = y_var, latentVars = c(u_var, "mu"),
                     data_ind,
                     s_v, s_p, t_p, h_v, h_p, e_v, d, sig, sigsq)
  mod <- mxModel("mod", mod_dir, mod_ind, 
                 mxFitFunctionMultigroup(c("dir.fitfunction","ind.fitfunction")),
                 b, d, inc, dir_est, ind_est, 
                 mxCI(c('dir_est','ind_est','inc'), interval = conf))
  
  res <- mxRun(mod, intervals = T)
  output <- res$output$confidenceIntervals[, c(2,1,3)]
  rownames(output) <- c("Direct Estimate", "Indirect Estimate", "Difference")
  colnames(output) <- c("Estimate", "Lower Bound", "Upper Bound")
  if(printout){
    cat("Evidence-splitting model for", comp[2], "compared to", comp[1], ":\n")
    print(output)
  }
  return(res)
}

#### Dong2013 LABA-ICS compared to LABA ####

data(Dong2013, package = "netmeta")

Dong_new <- as_tibble(Dong2013) %>% 
  group_by(id) %>% 
  summarize(have_zero = (prod(death) == 0) | (prod(randomized - death) ==0), 
            all_zero = (sum(death) == 0) | (sum(randomized - death) ==0)) %>% 
  filter(!all_zero) %>%
  inner_join(Dong2013) %>%
  mutate(y = log((death + 0.5 * have_zero) / (randomized - death + 0.5 * have_zero)),
         v = 1/(death + 0.5 * have_zero) + 1/(randomized - death + 0.5 * have_zero))

Dong_res <- nma_incon(y = "y", v = "v", trt = "treatment", s = "id", 
                      comp = c("LABA", "LABA-ICS"), data = Dong_new)

#### Linde2015 Hypericum compared to Placebo ####

data(Linde2015, package = "netmeta")

Linde_new <- as.tibble(Linde2015) %>%
  select(id, starts_with("treatment"), starts_with("resp"), starts_with("n")) %>%
  pivot_longer(treatment1:n3, names_to = ".value", names_pattern = "(.*).") %>%
  filter(!is.na(resp))

Linde_new <- Linde_new %>% 
  group_by(id) %>% 
  summarize(have_zero = (prod(resp) == 0) | (prod(resp - n) ==0), 
            all_zero = (sum(resp) == 0) | (sum(resp - n) ==0)) %>% 
  filter(!all_zero) %>%
  inner_join(Linde_new) %>%
  mutate(y = log((resp + 0.5 * have_zero) / (n - resp + 0.5 * have_zero)),
         v = 1/(resp + 0.5 * have_zero) + 1/(n - resp + 0.5 * have_zero))

Linde_res <- nma_incon(y = "y", v = "v", trt = "treatment", s = "id", 
                      comp = c("Placebo", "Hypericum"), data = Linde_new)

### Pagliaro1992 (the data used in Shih & Tu, 2021) ###
### Beta-blocker compared to Sclerotherapy          ###

data(Pagliaro1992 , package = "meta")
Pagliaro_new <- as_tibble(Pagliaro1992) %>% 
  group_by(id) %>% 
  slice(1) %>% 
  ungroup() %>%
  select(id, bleed.plac, n.plac) %>% 
  rename(bleed.exp = bleed.plac, n.exp = n.plac) %>%
  mutate(treat.exp = "Placebo") %>%
  bind_rows(as_tibble(Pagliaro1992) %>% 
              select(id, treat.exp, bleed.exp, n.exp))

Pagliaro_new <- Pagliaro_new %>% 
  group_by(id) %>% 
  summarize(have_zero = (prod(bleed.exp) == 0) | (prod(n.exp - bleed.exp) ==0), 
            all_zero = (sum(bleed.exp) == 0) | (sum(n.exp - bleed.exp) ==0)) %>% 
  filter(!all_zero) %>%
  inner_join(Pagliaro_new) %>%
  mutate(y = log((bleed.exp + 0.5 * have_zero) / (n.exp - bleed.exp + 0.5 * have_zero)),
         v = 1/(bleed.exp + 0.5 * have_zero) + 1/(n.exp - bleed.exp + 0.5 * have_zero))

Pagliaro_res <- nma_incon(y = "y", v = "v", trt = "treat.exp", s = "id", 
                          comp = c("Sclerotherapy", "Beta-blocker"), data = Pagliaro_new)
