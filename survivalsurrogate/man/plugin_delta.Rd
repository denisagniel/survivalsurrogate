\name{plugin_delta}
\alias{plugin_delta}
\title{
Plug-in Estimator of Causal Effect Based on EIF
}
\description{
Estimates a causal parameter using a plug-in estimator that combines machine learning-based nuisance function estimation and efficient influence function (EIF) correction.
}
\usage{
plugin_delta(data, folds, id, x, g, a = NULL, y, s, binary_lrnr = NULL, 
             cont_lrnr = NULL, e = NULL, gamma1 = NULL, gamma0 = NULL, 
             mu1 = NULL, mu0 = NULL, Q1 = NULL, Q0 = NULL, 
             Qstar1 = NULL, Qstar0 = NULL, truncate_e = 1e-12, 
             verbose = FALSE)
}
\arguments{
  \item{data}{A \code{data.table} or \code{data.frame} containing the observed data.}
  \item{folds}{A vector or list indicating cross-validation folds, used for sample splitting.}
  \item{id}{The name of the ID variable (as a string).}
  \item{x}{Vector of covariate column names used in nuisance estimation.}
  \item{g}{The name of the treatment indicator variable (as a string).}
  \item{a}{(Optional) Names of indicators for availability or adherence, if relevant.}
  \item{y}{Names of the outcome variables (can be longitudinal, e.g., \code{Y_1}, \code{Y_2}, ...).}
  \item{s}{Names of survival indicator variables, aligned with \code{y}.}
  \item{binary_lrnr}{Learner object used for binary outcome nuisance functions (e.g., treatment, censoring, hazard).}
  \item{cont_lrnr}{Learner object used for continuous outcome nuisance functions.}
  \item{e}{(Optional) Name of estimated propensity score variable. If NULL, it is estimated from the data.}
  \item{gamma1}{(Optional) Names of estimated censoring probabilities under treatment. If NULL, they are estimated.}
  \item{gamma0}{(Optional) Names of estimated censoring probabilities under control. If NULL, they are estimated.}
  \item{mu1}{(Optional) Names of estimated hazards under treatment. If NULL, they are estimated.}
  \item{mu0}{(Optional) Names of estimated hazards under control. If NULL, they are estimated.}
  \item{Q1}{(Optional) Names of estimated conditional means under treatment. If NULL, they are estimated.}
  \item{Q0}{(Optional) Names of estimated conditional means under control. If NULL, they are estimated.}
  \item{Qstar1}{(Optional) Placeholder for TMLE-style fluctuation targets (currently unused).}
  \item{Qstar0}{(Optional) Placeholder for TMLE-style fluctuation targets (currently unused).}
  \item{truncate_e}{Threshold for truncating small propensity scores to stabilize estimation.}
  \item{verbose}{Logical; if TRUE, prints messages during estimation.}
}
\details{
This function estimates a target causal effect using a plug-in approach based on nuisance functions: hazard, censoring, and outcome models. It also constructs the efficient influence function (EIF) and returns an estimate, its standard error, and the full EIF sample.
}
\value{
A \code{data.table} with the following components:
  \item{plugin_est}{Estimated causal effect.}
  \item{plugin_se}{Estimated standard error based on the EIF.}
  \item{if_data}{A nested \code{data.table} containing the EIF values for each observation.}
}
\references{
Agniel D and Parast L (2025+). "Robust Evaluation of Longitudinal Surrogate Markers with Censored Data." Journal of the Royal Statistical Society: Series B, In press.
}
