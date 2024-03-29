% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nmab.R
\name{nmab_gi}
\alias{nmab_gi}
\title{Calculate the Gittins index for a single arm (normal rewards)}
\usage{
nmab_gi(Sigma, n, gamma, tau, N, xi, delta, tol = 5e-04, lb = NA, ub = NA)
}
\arguments{
\item{Sigma}{Numeric. Value of Sigma for the arm (Bayesian reward).}

\item{n}{Numeric > 0. Value of n for the arm (Bayesian number of observation).}

\item{gamma}{Numeric in (0, 1). Reward discount factor.}

\item{tau}{Numeric > 0. Observation precision.}

\item{N}{Integer >= 2. Time horizon used.}

\item{xi}{Numeric > 0. Value of xi (entent of dynamic program state space).}

\item{delta}{Numeric > 0. Value of delta (fineness of discretisation in the dynamic program).}

\item{tol}{Numeric > 0. Absolute accuracy required.}

\item{lb}{Optional lower bound for GI.}

\item{ub}{Optional upper bound for GI.}
}
\value{
A vector of GI values.
}
\description{
The problem state is given by the \code{Sigma}, \code{n}, \code{gamma}, and \code{tau} arguments. The
remaining arguments affect how the calculation is done and are chosen based on accuracy
and speed requirements (see details).
}
\details{
The problem has an infinite continuous state space, but the calculation can only be done for a
finite discrete area. \code{N} and \code{xi} control the finite extent of the state space used (one in
each direction). These have a diminishing effect on accuracy and just need to be sufficiently
large, beyond which there will be minimal improvement. Each parameter has a linear effect on the state
space size and therefore algorithm speed. Larger \code{gamma} or smaller \code{tau} requires
a larger \code{N}, but depends only on the least favourable of the two. Often \code{N} can be smaller
than for similar BMAB problems. The required \code{xi} is more robust to changes in the problem setting
as adjusts to \code{tau} somewhat. Therefore, standard values can be used across problems with less
disadvantage. \code{xi} is based on standard deviations so 3 or 4 will often be large enough. Using \code{N} or
\code{xi} values that are too low gives an underestimate of the Gittins index.

The discretisation parameter \code{delta} is the main limiter to accuracy as it has a non-linear
effect on computation time. Accuracy increases as \code{delta} gets smaller, with inaccuracies leading
to an overestimate of the Gittins index.

The \code{lb} and \code{ub} arguments can be used to provide a starting interval for calibration if desired.
However, for normal use this is not needed as they will be calculated internally if not supplied.
So the initial interval is determined as follows:
\itemize{
\item For lower bound, use \code{lb} if supplied else use KGI.
\item For upper bound, use \code{ub} if supplied else use GI+.
}
}
\seealso{
For a link to the accompanying paper see \link{gittins-package}.
}
