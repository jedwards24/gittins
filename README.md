# Gittins

This is an R package to calculate Gittins indices for the multi-armed bandit problem.

## Description

This project contains functions written in R to calculate Gittins indices for the Bayesian multi-armed bandit problem with Bernoulli or Normal rewards. More information on the methodology can be found in this [paper](https://arxiv.org/abs/1909.05075) and in my [thesis](http://eprints.lancs.ac.uk/84589/) (in Chapter 5).

## Installation

Install the gittins package from github with:

``` r
# install.packages("devtools")
devtools::install_github("jedwards24/gittins")
```

## Usage

There are two main groups of functions, those for the multi-armed bandit with bernoulli rewards (begin with `bmab_`) and those for the same problem with normal rewards (begin with `nmab_`). The functions `bmab_gi()` and `nmab_gi()` will calculate the gittins index of a single arm state, however, it is more common to calculate indices for the multiple states that will be encountered in a multi-armed bandit problem.

### Bernoulli Rewards

Give the starting state of the problem (using two parameters, `alpha` and `beta`) and the number of actions that will be taken. These define the possible states the arm could take and therefore the states for which Gittins indices are needed, which are returned by the following function as a triangular matrix:

``` r
bmab_gi_multiple(alpha_start = 1, beta_start = 1, gamma = 0.9, N = 80, num_actions = 20, tol = 5e-5)
```

The other arguments used are:  the discount factor `gamma` of the problem; the desired accuracy `tol`; and the horizon `N` of the dynamic programme. The first of these come from the problem while the last two affect the calculation accuracy and speed. The alternative parameterisation of `Sigma` and `n` can be used by giving `Sigma_start` and `n_start` arguments instead of `alpha_start` and `beta_start`.

### Normal Rewards

For this problem, Gittins indices are only needed for a single vector of states defined by the parameter `n`. These are given by `n_range` (which must be ascending). The following function calculates indices assuming our mean belief in reward of the arms is zero (indices for different values can be found by a simple transform):

``` r
nmab_gi_multiple(n_range = 1 : 20, gamma = 0.9, tau = 1, tol = 5e-5, N = 30, xi = 3, delta = 0.02)
```

The arguments `gamma` and `tau` are problem settings while the remaining arguments affect the solution accuracy. Calculations for the NMAB can be much slower than for the BMAB especially if high accuracy is required.
