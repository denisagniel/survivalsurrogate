\name{estimate_R}
\alias{estimate_R}

\title{
Estimates the proportion of the treatment effect explained
}
\description{
Estimates the proportion of the treatment effect on the censored primary outcome that is explained by the longitudinal surrogate marker, given the estimated influence functions for the overall treatment effect (delta) and the residual treatment effect (delta_s).
}
\usage{
estimate_R(delta_if, delta_s_if, delta = NULL, delta_s = NULL, se_type = "asymptotic",
n_boot = NULL, alpha = 0.05)
}
\arguments{
  \item{delta_if}{
Influence function estimates for delta.
}
  \item{delta_s_if}{
Influence function  estimates for delta_s.
}
  \item{delta}{
(Optional) Defaults to the mean of the delta influence function estimates.}
  \item{delta_s}{
(Optional) Defaults to the mean of the delta_s influence function estimates.
}
  \item{se_type}{
Type of standard error estimation, choices are "asymptotic" or "bootstrap".}
  \item{n_boot}{
(Optional unless se_type = "bootstrap") Number of boostrap samples.
}
  \item{alpha}{
(Optional) Alpha level used for confidence interval. If not provided, this is set to 0.05.
}
}

\value{
A dataframe with the following components:
  \item{estimand}{Indicates what estimand is estimated, including the treatment effect (delta), residual treatment effect (delta.s), and proportion of treatment effect explained (R) by the longitudinal surrogate marker.}
  \item{estimate}{Estimates of the treatment effect (delta), residual treatment effect (delta.s), and proportion of treatment effect explained (R) by the longitudinal surrogate marker.}
  \item{se}{Estimated standard errors.}
  \item{ci_l}{Lower bound of the confidence intervals.}
   \item{ci_h}{Upper bound of the confidence intervals.}
}
\references{
Agniel D and Parast L (2025). "Robust Evaluation of Longitudinal Surrogate Markers with Censored Data." Journal of the Royal Statistical Society: Series B; \doi{10.1093/jrsssb/qkae119}.
}

\examples{
data(exampledata)
names(exampledata)
library(glue)
library(rpart)
library(mlr3)
library(dplyr)
library(mlr3learners)

tt <- 5
t0 <- 4
yvars <- paste0('Y_', 0:tt)
lrnc <- glue('regr.rpart')
lrnb <- glue('classif.log_reg')

\donttest{
p_deltahat <- plugin_delta(
    data = exampledata,
    folds = 'ff',
    id = 'ID',
    x = 'X_0',
    g = 'G_0',
    a = paste0('A_', 0:tt),
    y = yvars,
    s = paste0('S_', 0:t0),
    binary_lrnr = lrn(lrnb, predict_type = 'prob'),
    cont_lrnr = lrn(lrnc),
    truncate_e = 0.005,
    verbose = FALSE
  )

p_deltahat_s <- plugin_delta_s(
    data = exampledata,
    folds = 'ff',
    id = 'ID',
    x = 'X_0',
    g = 'G_0',
    a = paste0('A_', 0:tt),
    y = yvars,
    s = paste0('S_', 0:t0),
    binary_lrnr = lrn(lrnb, predict_type = 'prob'),
    cont_lrnr = lrn(lrnc),
    t0=t0,
    truncate_pi = 0.005,
    truncate_e = 0.005,
    verbose = FALSE
  )


estimate_R(p_deltahat$if_data[[1]]$eif,
p_deltahat_s$if_data[[1]]$eif)
}
}
