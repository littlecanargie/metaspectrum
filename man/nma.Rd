% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nma.R
\name{nma}
\alias{nma}
\title{Network meta-analysis under structural equation modelling}
\usage{
nma(
  s = "s",
  y = "y",
  v = "v",
  trt = "trt",
  trt_list = NULL,
  baseline = FALSE,
  trt_het = FALSE,
  uwls = FALSE,
  data = NULL,
  printout = TRUE,
  large_var = 1e+05,
  conf = 0.95
)
}
\arguments{
\item{s}{Column name for study (default "s")}

\item{y}{Column name for reported outcome (default "y")}

\item{v}{Column name for outcome variance (default "v")}

\item{trt}{Column name for treatment (default "trt"); s, y, v, trt
can be input as vectors directly if data is set to NULL}

\item{trt_list}{Vector of treatment in the desired order for the
output}

\item{baseline}{Set the study effect as fixed or random (default
FALSE)}

\item{trt_het}{Distinct degree of heterogeneity for each treatment
(default FALSE)}

\item{uwls}{Use unweighted least squares (default FALSE)}

\item{data}{Input data frame (long form, i.e. one study-treatment
pair per row)}

\item{printout}{Boolean for printing out the results (default TRUE)}

\item{large_var}{Large variance for incidental parameter problem
(default 1e5)}

\item{conf}{confidence level of the confidence intervals (default
.95); confidence intervals are constructed via likelihood-based
approach}
}
\value{
List ...
}
\description{
Network meta-analysis under structural equation modelling
}
\details{
Network meta-analysis under structural equation modelling

!! Currently only support continuous outcomes !!

!! Currently does not check if the network is connected!!
}
\examples{
library("dplyr")
library("tidyr")

#### Dong2013 ####

data(Dong2013, package = "netmeta")

Dong_new <- as_tibble(Dong2013) \%>\% 
  group_by(id) \%>\% 
  summarize(have_zero = (prod(death) == 0) | (prod(randomized - death) ==0), 
            all_zero = (sum(death) == 0) | (sum(randomized - death) ==0)) \%>\% 
  filter(!all_zero) \%>\%
  inner_join(Dong2013) \%>\%
  mutate(y = log((death + 0.5 * have_zero) / (randomized - death + 0.5 * have_zero)),
         v = 1/(death + 0.5 * have_zero) + 1/(randomized - death + 0.5 * have_zero))

Dong_res <- nma(y = "y", v = "v", trt = "treatment", s = "id", data = Dong_new)

#### Linde2015 ####

data(Linde2015, package = "netmeta")

Linde_new <- as_tibble(Linde2015) \%>\%
  select(id, starts_with("treatment"), starts_with("resp"), starts_with("n")) \%>\%
  pivot_longer(treatment1:n3, names_to = ".value", names_pattern = "(.*).") \%>\%
  filter(!is.na(resp))

Linde_new <- Linde_new \%>\% 
  group_by(id) \%>\% 
  summarize(have_zero = (prod(resp) == 0) | (prod(resp - n) ==0), 
            all_zero = (sum(resp) == 0) | (sum(resp - n) ==0)) \%>\% 
  filter(!all_zero) \%>\%
  inner_join(Linde_new) \%>\%
  mutate(y = log((resp + 0.5 * have_zero) / (n - resp + 0.5 * have_zero)),
         v = 1/(resp + 0.5 * have_zero) + 1/(n - resp + 0.5 * have_zero))

Linde_res <- nma(y = "y", v = "v", trt = "treatment", s = "id", data = Linde_new)

### Pagliaro1992 (the data used in Shih & Tu, 2019) ###

data(Pagliaro1992 , package = "meta")
Pagliaro_new <- as_tibble(Pagliaro1992) \%>\% 
  group_by(id) \%>\% 
  slice(1) \%>\% 
  ungroup() \%>\%
  select(id, bleed.plac, n.plac) \%>\% 
  rename(bleed.exp = bleed.plac, n.exp = n.plac) \%>\%
  mutate(treat.exp = "Placebo") \%>\%
  bind_rows(as_tibble(Pagliaro1992) \%>\% 
              select(id, treat.exp, bleed.exp, n.exp))

Pagliaro_new <- Pagliaro_new \%>\% 
  group_by(id) \%>\% 
  summarize(have_zero = (prod(bleed.exp) == 0) | (prod(n.exp - bleed.exp) ==0), 
            all_zero = (sum(bleed.exp) == 0) | (sum(n.exp - bleed.exp) ==0)) \%>\% 
  filter(!all_zero) \%>\%
  inner_join(Pagliaro_new) \%>\%
  mutate(y = log((bleed.exp + 0.5 * have_zero) / (n.exp - bleed.exp + 0.5 * have_zero)),
         v = 1/(bleed.exp + 0.5 * have_zero) + 1/(n.exp - bleed.exp + 0.5 * have_zero))

Pagliaro_res <- nma(y = "y", v = "v", trt = "treat.exp", s = "id", 
                    trt_list = c("Placebo", "Sclerotherapy", "Beta-blocker"),
                    data = Pagliaro_new)

Pagliaro_het <- nma(y = "y", v = "v", trt = "treat.exp", s = "id", 
                    trt_list = c("Placebo", "Sclerotherapy", "Beta-blocker"),
                    trt_het = TRUE, data = Pagliaro_new)

Pagliaro_uwls <- nma(y = "y", v = "v", trt = "treat.exp", s = "id", 
                     trt_list = c("Placebo", "Sclerotherapy", "Beta-blocker"),
                     uwls= TRUE, data = Pagliaro_new)

Pagliaro_uwls_het <- nma(y = "y", v = "v", trt = "treat.exp", s = "id", 
                         trt_list = c("Placebo", "Sclerotherapy", "Beta-blocker"),
                         trt_het = TRUE, uwls= TRUE, data = Pagliaro_new)

}
\references{
Shih ...
Title.
\emph{Journal},
\bold{x}, p1--p2
}
\seealso{
\code{\link{nma_incon}}
}
\author{
Ming-Chieh Shih \email{mcshih@mx.nthu.edu.tw}
}
