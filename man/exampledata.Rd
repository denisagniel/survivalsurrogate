\name{exampledata}
\alias{exampledata}
\docType{data}
\title{
Example dataset
}
\description{
Example dataset with longitudinal surrogate marker and censored primary outcome}
\usage{data("exampledata")}
\format{
  A data frame with 1000 observations on the following 21 variables. Here, the landmark time is 4 and the final time point is 5. The surrogate is measured at the following time points: 0,1,2,3,4.
  \describe{
    \item{\code{ID}}{Unique individual ID}
    \item{\code{X_0}}{Baseline covariate vector.}
    \item{\code{G_0}}{Treatment indicator vector.}
    \item{\code{Y_0}}{Primary outcome variable at time 0.}
    \item{\code{S_0}}{Surrogate marker value at time 0.}
    \item{\code{Y_1}}{Primary outcome variable at time 1.}
    \item{\code{S_1}}{Surrogate marker value at time 1.}
    \item{\code{Y_2}}{Primary outcome variable at time 2.}
    \item{\code{S_2}}{Surrogate marker value at time 2.}
    \item{\code{Y_3}}{Primary outcome variable at time 3.}
    \item{\code{S_3}}{Surrogate marker value at time 3.}
    \item{\code{Y_4}}{Primary outcome variable at time 4.}
    \item{\code{S_4}}{Surrogate marker value at time 4.}
    \item{\code{Y_5}}{Primary outcome variable at time 5.}
    \item{\code{ff}}{Fold number.}
    \item{\code{A_0}}{Observation variable at time 0.}
    \item{\code{A_1}}{Observation variable at time 1.}
    \item{\code{A_2}}{Observation variable at time 2.}
    \item{\code{A_3}}{Observation variable at time 3.}
    \item{\code{A_4}}{Observation variable at time 4.}
    \item{\code{A_5}}{Observation variable at time 5.}
  }
}
