% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmab.R
\name{bmab_giplus}
\alias{bmab_giplus}
\title{Calculate the GI+ index for a single arm (Bernoulli rewards)}
\usage{
bmab_giplus(Sigma, n, gamma, tol = 5e-04, upper = FALSE)
}
\arguments{
\item{Sigma}{Numeric > 0. Value of Sigma for the arm.}

\item{n}{Numeric > Sigma > 0. Value of n for the arm.}

\item{gamma}{Numeric in (0, 1). Reward discount factor.}

\item{tol}{Numeric > 0. Absolute accuracy required.}

\item{upper}{if TRUE, the upper end of the interval is returned, otherwise the midpoint.}
}
\value{
A GI+ index value.
}
\description{
The GI+ index is an upper bound for the Gittins index.
}
