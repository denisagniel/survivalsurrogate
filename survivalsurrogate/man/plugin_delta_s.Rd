\name{plugin_delta_s}
\alias{plugin_delta_s}
\title{
Plugin Estimator of the Surrogate Causal Effect at a Fixed Time Point
}
\description{
Computes a plugin estimator of the average causal effect of treatment on a primary outcome, conditionally mediated through a surrogate marker, at a fixed time point. This estimator uses flexible machine learning models for nuisance parameters and supports both asymptotic and bootstrap standard errors.
}
\usage{
plugin_delta_s(data, folds, id, x, g, a = NULL, y, s, binary_lrnr = NULL, cont_lrnr = NULL, t0 = length(s), e = NULL, gamma1 = NULL, gamma0 = NULL, mu1 = NULL, mu0 = NULL, pi = NULL, pistar = NULL, Q1 = NULL, Q0 = NULL, Qstar1 = NULL, Qstar0 = NULL, truncate_e = 1e-12, truncate_pi = 1e-12, se_type = "asymptotic", n_boot = NULL, alpha = 0.05, verbose = FALSE, retain_data = FALSE)
}
\arguments{
  \item{data}{A \code{data.table} or \code{data.frame} with one row per observation.}
  \item{folds}{A vector of fold assignments used for cross-fitting.}
  \item{id}{The name of the column identifying individual observations.}
  \item{x}{A character vector of baseline covariate column names.}
  \item{g}{The name of the treatment assignment column (0 or 1).}
  \item{a}{Optional; vector of censoring indicator column names at each time point. If \code{NULL}, assumes no censoring.}
  \item{y}{Vector of primary outcome indicator column names at each time point.}
  \item{s}{Vector of surrogate marker indicator column names at each time point.}
  \item{binary_lrnr}{A binary outcome learner (e.g., a classification model) used to estimate binary nuisance parameters such as hazards and censoring probabilities.}
  \item{cont_lrnr}{A continuous outcome learner (e.g., regression model) used to estimate conditional expectations.}
  \item{t0}{The time index at which the causal effect is evaluated. Defaults to the last time point.}
  \item{e}{Optional; the name of the column containing estimated treatment propensity scores. If \code{NULL}, they are estimated.}
  \item{gamma1}{Optional; vector of column names for estimated censoring probabilities under treatment. If \code{NULL}, they are estimated or assumed to be 1.}
  \item{gamma0}{Optional; vector of column names for estimated censoring probabilities under control. If \code{NULL}, they are estimated or assumed to be 1.}
  \item{mu1}{Optional; vector of column names for hazards under treatment. If \code{NULL}, they are estimated.}
  \item{mu0}{Optional; vector of column names for hazards under control. If \code{NULL}, they are estimated.}
  \item{pi}{Optional; vector of column names for estimated probabilities of surrogate marker under control. If \code{NULL}, they are estimated.}
  \item{pistar}{Optional; vector of column names for estimated probabilities of surrogate marker under treatment. If \code{NULL}, they are estimated.}
  \item{Q1}{Optional; vector of column names for conditional expectations of outcome under treatment. If \code{NULL}, they are estimated.}
  \item{Q0}{Optional; vector of column names for conditional expectations of outcome under control. If \code{NULL}, they are estimated.}
  \item{Qstar1}{Optional; vector of column names for pseudo-outcome expectations under treatment. If \code{NULL}, they are estimated.}
  \item{Qstar0}{Optional; vector of column names for pseudo-outcome expectations under control. If \code{NULL}, they are estimated.}
  \item{truncate_e}{Truncation threshold for propensity scores to avoid division by zero. Default is \code{1e-12}.}
  \item{truncate_pi}{Truncation threshold for surrogate probabilities to avoid instability. Default is \code{1e-12}.}
  \item{se_type}{Method to estimate standard errors: either \code{"asymptotic"} or \code{"bootstrap"}.}
  \item{n_boot}{Number of bootstrap replicates if \code{se_type = "bootstrap"}. Required if using bootstrap.}
  \item{alpha}{Significance level for confidence intervals. Default is \code{0.05}.}
  \item{verbose}{Logical; if \code{TRUE}, prints progress messages.}
  \item{retain_data}{Logical; if \code{TRUE}, retains full dataset with estimated components in the output.}
}
\details{
This function computes a plug-in estimator for the surrogate causal effect, i.e., the average causal effect of treatment on the primary outcome that is mediated by a surrogate marker, at a specific time point. The function optionally estimates nuisance parameters (e.g., hazards, expectations, surrogate probabilities) using cross-fitting and machine learning, and returns estimates with standard errors and confidence intervals.
}
\value{
A \code{data.frame} with the following components:
  \item{plugin_est}{The plug-in point estimate of the surrogate causal effect at time \code{t0}.}
  \item{plugin_se}{Standard error of the estimate, based on influence function or bootstrap.}
  \item{ci_l}{Lower bound of the (1 - \code{alpha}) confidence interval.}
  \item{ci_h}{Upper bound of the (1 - \code{alpha}) confidence interval.}
  \item{if_data}{A list containing the dataset with the efficient influence function if \code{retain_data = TRUE}.}
}
\references{
Agniel D and Parast L (2025+). "Robust Evaluation of Longitudinal Surrogate Markers with Censored Data." Journal of the Royal Statistical Society: Series B, In press.
}
