#' @description
#' This is an R package to calculate Gittins indices for the multi-armed bandit problem.
#'
#' This project contains functions written in R to calculate Gittins indices for
#' the Bayesian multi-armed bandit problem with Bernoulli or Normal rewards. Detailed information
#' on the methodology as well as the effect of different computation parameters on speed and
#' accuracy can be found in the paper at <https://arxiv.org/abs/1909.05075>
#' and in my thesis at <http://eprints.lancs.ac.uk/84589/> (in Chapter 5).
#'
#' I welcome any feedback or suggestions (see link to Github page below). If the package is useful
#'  in your work then can you please reference the package and/or the paper.
#'
#' @author James Edwards (also maintainer)
#' @seealso
#' * The main package Github page is <https://github.com/jedwards24/gittins>.
#' * Report bugs at <https://github.com/jedwards24/gittins/issues>.
#'
#' @importFrom stats dbeta dnorm integrate pbeta pnorm na.omit
#' @importFrom utils setTxtProgressBar txtProgressBar
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
