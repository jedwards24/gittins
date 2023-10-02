# Gittins

This is an R package to calculate Gittins indices for the multi-armed bandit problem.

## Description

This project contains functions written in R to calculate Gittins indices for the Bayesian multi-armed bandit problem with Bernoulli or Normal rewards. Detailed information on the methodology as well as the effect of different computation parameters on speed and accuracy can be found in this [paper](https://arxiv.org/abs/1909.05075) and in my [thesis](http://eprints.lancs.ac.uk/84589/) (in Chapter 5).

## Installation

Install the gittins package from github with:

``` r
# install.packages("devtools")
devtools::install_github("jedwards24/gittins")
```

## Usage

There are two main groups of functions, those for the multi-armed bandit with bernoulli rewards (begin with `bmab_`) and those for the same problem with normal rewards (begin with `nmab_`). 

The functions `bmab_gi()` and `nmab_gi()` will calculate the gittins index of a single arm state, while `bmab_gi_multiple()` and `nmab_gi_multiple()` will calculate indices for the multiple states that will be encountered in a multi-armed bandit problem.

The package also contains functions to calculate other indices, the KGI and GI+, which provide lower and upper bounds 
for the Gittins index.

### Bernoulli Rewards

The following calculates the Gittins index for an arm which has a Bayesian state (prior plus observed) of 
two observations and one success. 

``` r
bmab_gi(1, 2, gamma = 0.9, N = 80, tol = 5e-5)
```

The other arguments used are:  the reward discount factor `gamma` of the problem, the desired accuracy `tol`, and the horizon `N` of the dynamic program. The first of these is a feature of the problem while the last two affect the calculation accuracy and speed. 

The alternative function `bmab_gi_ab()` gives the arm state in alpha/beta values (success/failure counts). 

For multiple states, give the starting state of the problem (here using alpha/beta terminology) and the number of actions that will be taken. These define the possible states the arm could take and therefore the states for which Gittins indices are needed. The following function gives these results in a data frame:

``` r
bmab_gi_multiple(alpha_start = 1, beta_start = 1, num_actions = 10, gamma = 0.9, N = 80, tol = 5e-5)
```

### Normal Rewards

This is similar to the Bernoulli case but with extra parameters. For a a single state:

```r
nmab_gi(0, 1, gamma = 0.9, tau = 1, N = 30, xi = 3, delta = 0.02, tol = 5e-4)
```

The arguments `gamma` and `tau` are problem setting parameters, while the remaining arguments affect 
the solution accuracy and speed. Calculations for the NMAB can be much slower than for the BMAB 
especially if high accuracy is required.

For multiple states, Gittins indices are only needed for a single vector of states given by the parameter `n` since indices for different mean beliefs can then be can be found by a simple transform. 

``` r
nmab_gi_multiple(n_range = 1:10, gamma = 0.9, tau = 1, N = 30, xi = 3, delta = 0.02, tol = 5e-4)
```
