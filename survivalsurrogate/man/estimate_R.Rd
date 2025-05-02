\name{estimate_R}
\alias{estimate_R}

\title{
Estimates the proportion of the treatment effect explained
}
\description{
Estimates the proportion of the treatment effect on the censored primary outcome that is explained by the longitudinal surrogate marker
}
\usage{
estimate_R(delta_if, delta_s_if, delta = NULL, delta_s = NULL, se_type = "asymptotic", n_boot = NULL, alpha = 0.05)
}
\arguments{
  \item{delta_if}{
delta influence function
}
  \item{delta_s_if}{
delta.s influence function
}
  \item{delta}{
optional; defaults to the mean of the delta influence function}
  \item{delta_s}{
optional; defaults to the mean of the delta.s influence function
}
  \item{se_type}{
choices are "asymptotic" or "bootstrap"}
  \item{n_boot}{
number of boostrap samples, required if se_type = "bootstrap"
}
  \item{alpha}{
alpha level used for confidence interval, defaults to 95\% confidence interval
}
}

\value{
returns a tibble with the estimate, standard error (se), and confidence interval limits for delta, delta.s, and R (the proportion of treatment effect explained)
}
\references{
Agniel D and Parast L (2025+). "Robust Evaluation of Longitudinal Surrogate Markers with Censored Data." Journal of the Royal Statistical Society: Series B, In press.
}
