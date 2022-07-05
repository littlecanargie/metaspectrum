# =========================================================================== #
#        Network meta-analysis under structural equation modelling            #
# --------------------------------------------------------------------------- #
#             !! Currently only support continuous outcomes !!                #
#         !! Currently does not check if the network is connected!!           #
#                                                                             #
#  main function: nma(s, y, v, trt, comp, data, printout, large_var)          #
#                                                                             #
#  data: Input data frame (long form, i.e. one study-treatment pair per row)  #
#  s: Column name for study (default "s")                                     #
#  y: Column name for reported outcome (default "y")                          #
#  v: Column name for outcome variance (default "v")                          #
#  trt: Column name for treatment (default "trt")                             #
#    -- s, y, v, trt can be input as vectors directly if data is set to NULL  #
#  trt_list: Vector of treatment in the desired order for the output          #
#  baseline: Set the study effect as fixed or random (default F)              #
#  trt_het: Distinct degree of heterogeneity for each treatment (default F)   #
#  uwls: Use unweighted least squares (default F)                             #
#  printout: Boolean for printing out the results (default T)                 #
#  large_var: Large variance for incidental parameter problem (default 1e5)   #
#  conf: confidence level of the confidence intervals (default .95)           #
#    -- Confidence intervals are constructed via likelihood-based approach    #
#                                                                             #
# --------------------------------------------------------------------------- #
#  Author : Ming-Chieh Shih (Last updated 2022-07-05)                         #
# =========================================================================== #

require.i <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# netmeta and meta are called for real-world data 
# (Dong2013, Linde 2015, Pagliaro1992)
require.i(c("OpenMx", "tidyverse", "netmeta", "meta"))

nma <- function(s = "s", y = "y", v = "v", trt = "trt", trt_list = NULL,
                baseline = F, trt_het = F, uwls = F,
                data = NULL, printout = T, large_var = 1e5, conf = .95){
  if(is.null(data)) data <- sys.frame(sys.parent())
  
  # Obtain coding book for treatment recoding indices and number of treatments
  if(is.null(trt_list)) trt_list <- levels(factor(data[[trt]]))
  ntrt <- length(trt_list)
  
  # Pivot to wide table, define direct study with coi, and 
  # fill in placeholder variance for missing arms
  data <- data.frame(s = data[[s]], 
                     y = data[[y]], 
                     v = data[[v]], 
                     trt = as.numeric(factor(data[[trt]], levels = trt_list)))
  
  data <- data %>% 
    pivot_wider(id_cols = s, names_from = trt, 
                names_sep = "", values_from = c(y, v)) %>%
    mutate(across(starts_with('v'), ~replace_na(., 1)))
  
  # Variable name vectors for future use
  y_var <- paste0("y", 1:ntrt)
  u_var <- paste0("u", 1:ntrt)
  v_var <- paste0("data.v", 1:ntrt)
  e_var <- paste0("e", 1:ntrt)
  
  # Definition of parameter matrices
  # d (d1 - d'ntrt'): column vector of size ntrt for treatment effects 
  # sigsq (sigsq): half the heterogeneity variance
  # sqrtphimat (sqrtphi): square root of the uwls weight 
  d <- mxMatrix(type = "Full", nrow = ntrt, ncol = 1, free = c(F, rep(T, ntrt-1)),
                labels = c(paste0("d", 1:ntrt)), values = rep(0, ntrt), name = "dmat")
  
  sig_lb <- "sig"; if(trt_het) sig_lb <- paste0("sig", 1:ntrt)
  sig <- mxMatrix(type = "Full", nrow = ifelse(trt_het, ntrt, 1), ncol = 1, free = T, 
                  labels = sig_lb, name = "sigmat")
  sigsq <- mxAlgebra(sigmat^2, name = "sigsq")
  
  sqrtphi_lb <- "sqrtphi"; if(trt_het) sqrtphi_lb <- paste0("sqrtphi", 1:ntrt)
  sqrtphi <- mxMatrix(type = "Full", nrow = ifelse(trt_het, ntrt, 1), ncol = 1, free = T, 
                      labels = sqrtphi_lb, name = "sqrtphimat")
  
  # Model lines for study effects
  s_v <- mxPath(from = "mu", arrows = 2, free = baseline, 
                labels = "sig_bsq", values = ifelse(baseline, 1, large_var))
  s_p <- mxPath(from = "mu", to = y_var, free = F, values = 1)
  
  # Model lines for treatment effects
  t_p <- mxPath(from = "one", to = y_var, free = c(F, rep(T, ntrt-1)), 
                labels = paste0("d", 1:ntrt))
  
  # Model lines for heterogeneity effects
  h_trt <- 1; if(trt_het) h_trt <- 1:ntrt
  h_v <- mxPath(from = u_var, arrows = 2, free = F, 
                labels = paste0("sigsq[",h_trt,",1]"))
  h_p <- mxPath(from = u_var, to = y_var, free = F, values = 1)
  
  # Model lines for error variance
  e_v <- mxPath(from = e_var, arrows = 2, free = F, labels = v_var)
  if(uwls){
    e_p <- mxPath(from = e_var, to = y_var, free = T, labels = sqrtphi_lb)
  }else{
    e_p <- mxPath(from = e_var, to = y_var, free = F, values = 1)
  }
  
  # Extra matrices for output
  cmat <- matrix(0, nrow = ntrt*(ntrt-1)/2, ncol = ntrt)
  tmat <- matrix("", nrow = ntrt*(ntrt-1)/2, ncol = 2)
  nowrow <- 1
  for(i in 1:(ntrt-1)){
    for(j in (i+1):ntrt){
      cmat[nowrow, i] <- -1; tmat[nowrow, 1] <- trt_list[i]
      cmat[nowrow, j] <- 1; tmat[nowrow, 2] <- trt_list[j]
      nowrow <- nowrow + 1
    }
  }
  colnames(tmat) <- c("reference", "treatment")
  
  b <- mxMatrix(type = "Full", nrow = ntrt*(ntrt-1)/2, ncol = ntrt, free = F,
                values = c(cmat), name = "b")
  est <- mxAlgebra(b %*% dmat, name = "est")
  
  # Putting together the model
  if(uwls){
    mod <- mxModel("LGM", type="RAM", manifestVars = y_var, latentVars = c(u_var, e_var, "mu"),
                   mxData(data.frame(data), type = "raw"),
                   # parameter matrices
                   d, sqrtphi,
                   # study effect
                   s_v, s_p,
                   # treatment effect
                   t_p,
                   # error
                   e_v, e_p,
                   b, est, mxCI('est'))
  }else{
    mod <- mxModel("LGM", type="RAM", manifestVars = y_var, latentVars = c(u_var, e_var, "mu"),
                   mxData(data.frame(data), type = "raw"),
                   # parameter matrices
                   d, sig, sigsq,
                   # study effect
                   s_v, s_p,
                   # treatment effect
                   t_p,
                   # heterogeneity effect
                   h_v, h_p,
                   # error
                   e_v, e_p,
                   b, est, mxCI('est')) 
  }
  
  res <- mxRun(mod, intervals = T)
  output <- as.tibble(res$output$confidenceIntervals[, c(2,1,3)])
  output <- cbind(tmat, output)
  print(output)
  return(res)
}

#### Dong2013 ####

data(Dong2013, package = "netmeta")

Dong_new <- as_tibble(Dong2013) %>% 
  group_by(id) %>% 
  summarize(have_zero = (prod(death) == 0) | (prod(randomized - death) ==0), 
            all_zero = (sum(death) == 0) | (sum(randomized - death) ==0)) %>% 
  filter(!all_zero) %>%
  inner_join(Dong2013) %>%
  mutate(y = log((death + 0.5 * have_zero) / (randomized - death + 0.5 * have_zero)),
         v = 1/(death + 0.5 * have_zero) + 1/(randomized - death + 0.5 * have_zero))

Dong_res <- nma(y = "y", v = "v", trt = "treatment", s = "id", data = Dong_new)

#### Linde2015 ####

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

Linde_res <- nma(y = "y", v = "v", trt = "treatment", s = "id", data = Linde_new)

### Pagliaro1992 (the data used in Shih & Tu, 2019) ###

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

Pagliaro_res <- nma(y = "y", v = "v", trt = "treat.exp", s = "id", 
                    trt_list = c("Placebo", "Sclerotherapy", "Beta-blocker"),
                    data = Pagliaro_new)

Pagliaro_het <- nma(y = "y", v = "v", trt = "treat.exp", s = "id", 
                    trt_list = c("Placebo", "Sclerotherapy", "Beta-blocker"),
                    trt_het = T, data = Pagliaro_new)

Pagliaro_uwls <- nma(y = "y", v = "v", trt = "treat.exp", s = "id", 
                     trt_list = c("Placebo", "Sclerotherapy", "Beta-blocker"),
                     uwls= T, data = Pagliaro_new)

Pagliaro_uwls_het <- nma(y = "y", v = "v", trt = "treat.exp", s = "id", 
                         trt_list = c("Placebo", "Sclerotherapy", "Beta-blocker"),
                         trt_het = T, uwls= T, data = Pagliaro_new)
