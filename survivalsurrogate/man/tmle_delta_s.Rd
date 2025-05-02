\name{tmle_delta_s}
\alias{tmle_delta_s}
\title{
TMLE for Longitudinal Surrogate Evaluation with Delta Statistic
}
\description{
Estimates a robust treatment effect using targeted minimum loss-based estimation (TMLE) for evaluating longitudinal surrogate markers under potential censoring. The method computes a delta-type estimand based on the difference between counterfactual outcomes, using influence-function-based standard errors.
}
\usage{
tmle_delta_s(data, folds, id, x, g, a = NULL, y, s, binary_lrnr = NULL, cont_lrnr = NULL, t0 = length(s), e = NULL, gamma1 = NULL, gamma0 = NULL, pi = NULL, pistar = NULL, truncate_e = 1e-12, verbose = FALSE, retain_data = FALSE)
}
\arguments{
  \item{data}{
A \code{data.frame} containing all necessary variables for estimation.
}
  \item{folds}{
A vector or list defining cross-validation folds for nuisance estimation.
}
  \item{id}{
A string giving the name of the column with unique unit IDs (e.g., subject ID).
}
  \item{x}{
A character vector of covariate names to adjust for in nuisance estimation.
}
  \item{g}{
A string indicating the column name of the treatment indicator (typically 0/1).
}
  \item{a}{
(Optional) A character vector of time-varying treatment columns. If \code{NULL}, assumes static treatment.
}
  \item{y}{
A character vector of outcome variable names (possibly longitudinal).
}
  \item{s}{
A character vector of surrogate marker variable names (possibly longitudinal).
}
  \item{binary_lrnr}{
Learner object or specification used for estimating binary nuisance components (e.g., propensity scores, censoring).
}
  \item{cont_lrnr}{
Learner object used for estimating continuous-valued outcome regressions.
}
  \item{t0}{
Time index (default is the final time point) at which the delta estimand is evaluated.
}
  \item{e}{
(Optional) Column name or vector of propensity score estimates. If \code{NULL}, these will be estimated.
}
  \item{gamma1}{
(Optional) Vector of column names for censoring probabilities under treatment. Estimated if \code{NULL}.
}
  \item{gamma0}{
(Optional) Vector of column names for censoring probabilities under control. Estimated if \code{NULL}.
}
  \item{pi}{
(Optional) Vector of column names giving estimated probabilities of observed surrogate values under observed treatment. Estimated if \code{NULL}.
}
  \item{pistar}{
(Optional) Vector of column names giving estimated probabilities under the reference distribution. Estimated if \code{NULL}.
}
  \item{truncate_e}{
Numeric truncation level for propensity scores to avoid division by near-zero values. Default is \code{1e-12}.
}
  \item{verbose}{
Logical; if \code{TRUE}, print progress messages.
}
  \item{retain_data}{
Logical; if \code{TRUE}, return the full data and influence function values.
}
}
\value{
A \code{data.frame} with the following components:
  \item{tmle_est}{TMLE estimate of the treatment effect based on surrogate delta statistic.}
  \item{tmle_se}{Estimated standard error using the empirical standard deviation of the influence function.}
  \item{if_data}{(Optional) Full dataset joined with estimated influence function contributions.}
}
\references{
Agniel D and Parast L (2025+). "Robust Evaluation of Longitudinal Surrogate Markers with Censored Data." \emph{Journal of the Royal Statistical Society: Series B}, In press.
}
