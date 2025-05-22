\name{plugin_delta}
\alias{plugin_delta}
\title{
Treatment effect estimation using plug-in estimator
}
\description{
Estimates the treatment effect using the plug-in estimator.
}
\usage{
plugin_delta(data, folds, id, x, g, a = NULL, y, s, binary_lrnr = NULL, 
             cont_lrnr = NULL, e = NULL, gamma1 = NULL, gamma0 = NULL, 
             mu1 = NULL, mu0 = NULL, Q1 = NULL, Q0 = NULL, 
             truncate_e = 1e-12, verbose = FALSE)
}
\arguments{
   \item{data}{
A dataframe containing all necessary variables for estimation. The functions in this package require the dataframe to be in a specific form. Therefore, you will need to reformat your dataset as follows. At a minimum, you will need variables in your dataframe that indicate (1) the folds for crossfitting; (2) a unique observation identifier; (3) baseline covariates, if there are none you must have a variable with all values equal to 1 to provide to the argument x below; (4) a variable indicating treatment group which should be 1 for treatment and 0 for control; (5) a set of variables that contain the surrogate marker value at each time point up to and including the landmark time, denoted t0; (6) a set of variables that indicate observation status at each time point where the surrogate marker is measured, in addition to the final time which is denoted as t; (7)  a set of variables that indicate primary outcome status at each time point where the surrogate marker is measured, in addition to the final time which is denoted as t. See the exampledata for an example where t0=4 and t=5.
}
  \item{folds}{
A vector defining crossfitting folds for nuisance estimation.
}
  \item{id}{
A string giving the name of the column with unique unit IDs (e.g., individual ID).
}
  \item{x}{
A character vector of covariate names to adjust for in nuisance estimation. At a minimum this must have one covariate equal to 1 for all individuals.
}
  \item{g}{
A string indicating the column name of the treatment indicator with 1 for treatment and 0 for control.
}
 \item{a}{
(Optional) A character vector of the column names indicating observation status at each time point. Specifically, A_t = 1 if the individual is still under observation (i.e., uncensored) at time t, meaning their outcome status could in principle be measured at that time. A value of 0 means the individual was censored prior to t. If not provided, assumed to be 1 for all time points.
}
\item{y}{
A character vector of the column names indicating primary outcome status at each time point. Specifically, Y_t = 1 if the primary event (e.g., failure, relapse, death) has not yet occurred by time t, and the individual is still at risk. A value of 0 means the primary event occurred on or before time t.
}

  \item{s}{
A character vector of the column names indicating the surrogate marker values measured at each time point, where the value is NA if the individual is no longer observable at that time point. 
}
   \item{binary_lrnr}{
Learner object or specification used for estimating binary nuisance components (e.g., propensity scores, censoring).
}
  \item{cont_lrnr}{
Learner object used for estimating continuous-valued outcome regressions.
}
  \item{e}{
(Optional) Column name or vector of propensity score estimates. If not provided, these will be estimated.
}
  \item{gamma1}{
(Optional) Vector of column names for censoring probabilities under treatment. If not provided, these will be estimated.
}
  \item{gamma0}{
(Optional) Vector of column names for censoring probabilities under control. If not provided, these will be estimated.
}
\item{mu1}{(Optional) Names of estimated hazards under treatment. If not provided, these will be estimated.}
  \item{mu0}{(Optional) Names of estimated hazards under control. If not provided, these will be estimated.}
  \item{Q1}{(Optional) Names of estimated conditional means under treatment. If not provided, these will be estimated.}
  \item{Q0}{(Optional) Names of estimated conditional means under control.If not provided, these will be estimated.}
  \item{truncate_e}{
Numeric truncation level for propensity scores to avoid division by near-zero values. Default is \code{1e-12}.
}
  \item{verbose}{
Logical; if TRUE, print progress messages.
}
}

\value{
A dataframe with the following components:
  \item{plugin_est}{Plug-in estimate of the treatment effect}
  \item{plugin_se}{Estimated standard error of the plug-in estimator, based on the influence function.}
  \item{if_data}{A nested data frame containing the influence function contributions for each observation.}
}
\references{
Agniel D and Parast L (2025+). "Robust Evaluation of Longitudinal Surrogate Markers with Censored Data." Journal of the Royal Statistical Society: Series B, In press.
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
p_deltahat
}

