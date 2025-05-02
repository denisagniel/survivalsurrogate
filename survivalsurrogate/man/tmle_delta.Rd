\name{tmle_delta}
\alias{tmle_delta}
\title{
Targeted Maximum Likelihood Estimation with Censoring Adjustment
}
\description{
Estimates a causal effect using Targeted Maximum Likelihood Estimation (TMLE) while accounting for time-dependent outcomes and potential censoring under treatment and control. This implementation supports flexible machine learning-based nuisance parameter estimation.
}
\usage{
tmle_delta(data, folds, id, x, g, a = NULL, y, s, binary_lrnr = NULL, cont_lrnr = NULL, e = NULL, gamma1 = NULL, gamma0 = NULL, truncate_e = 1e-12, verbose = FALSE)
}
\arguments{
  \item{data}{
A \code{data.frame} containing all relevant variables for estimation, including covariates, treatment, censoring, and outcome variables.
}
  \item{folds}{
An object specifying the cross-validation folds to use for nuisance parameter estimation. Typically a list of fold assignments.
}
  \item{id}{
A character string giving the name of the column identifying individual units in the data.
}
  \item{x}{
A character vector specifying the baseline covariates to adjust for.
}
  \item{g}{
A character string giving the name of the treatment variable (binary indicator for treatment vs. control).
}
  \item{a}{
(Optional) A character vector of names for the time-varying treatment indicators. If \code{NULL}, assumes no censoring.
}
  \item{y}{
A character vector of names for the outcome variables at each time point.
}
  \item{s}{
A character vector of names for the indicators of being at risk (uncensored) at each time point.
}
  \item{binary_lrnr}{
A binary learner (e.g., SuperLearner object or algorithm name) used to estimate the propensity score and censoring mechanisms.
}
  \item{cont_lrnr}{
A continuous outcome learner used to estimate the conditional expectations of the outcomes.
}
  \item{e}{
(Optional) A character string giving the name of the column with known or externally estimated treatment propensity scores. If \code{NULL}, these are estimated from the data.
}
  \item{gamma1}{
(Optional) A character vector of names of columns giving the censoring probabilities under treatment. If \code{NULL}, these are estimated.
}
  \item{gamma0}{
(Optional) A character vector of names of columns giving the censoring probabilities under control. If \code{NULL}, these are estimated.
}
  \item{truncate_e}{
A numeric value indicating the truncation level for the propensity score, used to avoid division by near-zero values.
}
  \item{verbose}{
Logical; if \code{TRUE}, prints progress and estimation details.
}
}
\details{
This function estimates the average treatment effect on a time-indexed outcome using TMLE, incorporating estimates of the treatment assignment mechanism, censoring processes, and outcome regressions. It supports estimation under censoring by estimating or supplying the censoring probabilities under both treatment and control. The user can optionally supply known values for the propensity score and censoring probabilities.
}
\value{
A data frame with the following components:
  \item{tmle_est}{The estimated causal effect (average difference in potential outcomes).}
  \item{tmle_se}{The estimated standard error of the TMLE estimator, based on the influence function.}
  \item{if_data}{A nested data frame containing the influence function contributions for each observation.}
}
\references{
Agniel D and Parast L (2025+). "Robust Evaluation of Longitudinal Surrogate Markers with Censored Data." Journal of the Royal Statistical Society: Series B, In press.

}

