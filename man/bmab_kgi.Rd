% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmab.R
\name{bmab_kgi}
\alias{bmab_kgi}
\title{Calculate the knowledge gradient index for a single arm (Bernoulli rewards)}
\usage{
bmab_kgi(Sigma, n, gamma)
}
\arguments{
\item{Sigma}{Numeric > 0. Value of Sigma for the arm.}

\item{n}{Numeric > Sigma > 0. Value of n for the arm.}

\item{gamma}{Numeric in (0, 1). Reward discount factor.}
}
\value{
An index value or a vector of values.
}
\description{
The KGI is an lower bound for the Gittins index.
}
\details{
This is an exact closed form calculation. Arguments \code{Sigma} and \code{n} can be supplied as vectors in
which case a vector of index values will be returned.
}
