% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmab.R
\name{bmab_gi_value}
\alias{bmab_gi_value}
\title{Value of a one-armed bandit (Bernoulli rewards)}
\usage{
bmab_gi_value(lambda, Sigma, n, gamma, N)
}
\arguments{
\item{lambda}{Reward from the known arm.}

\item{Sigma}{Value of Sigma for the unknown arm.}

\item{n}{Value of n for the unknown arm.}

\item{gamma}{Numeric in (0, 1). Reward discount factor.}

\item{N}{Integer >= 2. Time horizon used in the calculation..}
}
\value{
Difference in value between safe and unknown arms.
}
\description{
Value of a one-armed bandit (Bernoulli rewards)
}
\keyword{internal}
