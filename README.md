[![Build Status](https://travis-ci.org/mschauer/PointProcessInference.jl.svg?branch=master)](https://travis-ci.org/mschauer/PointProcessInference.jl)

# PointProcessInference.jl
Fast and scalable non-parametric Bayesian inference for Poisson point processes

## Introduction

Poisson point processes are among the basic modelling tools in many areas. Their probabilistic properties are determined by their intensity function, the density *λ*.
This Julia package implements our non-parametric Bayesian approach to estimation of the intensity function *λ* for univariate Poisson point processes. For full details see our preprint

-  S. Gugushvili, M. Schauer, F. van der Meulen, and P. Spreij. Fast and scalable non-parametric Bayesian inference for Poisson point processes. __[arXiv:1804.03616 [stat.ME]](https://arxiv.org/abs/1804.03616)__, 2018.

Intuitively, a univariate Poisson point processes *X*, also called a non-homogeneous Poisson process, can be thought of as random
scattering of points in the time interval *[0,T]*, where the way the scattering occurs is determined by the intensity function *λ*.
An example is the ordinary Poisson process, for which the intensity *λ*  is constant.

## Methodology

We infer the intensity function *λ* in a non-parametric fashion. The function *λ* is a priori modelled as piecewise constant. This is even more natural, if the data have been already binned,
as is often the case in, e.g., astronomy. Thus, fix a positive integer *N* and a grid *b* of points `b[1] == 0`, `b[N] == T` on the interval *[0,T]*, for instance a uniform grid.
The intensity *λ* is then modelled as
`λ(x) = ψ[k]` for `b[k] <= x < b[k+1]`.

Now we postulate that a priori the coefficients *ψ* form a Gamma Markov chain (GMC). As explained in our preprint, this prior induces smoothing across the coefficients *ψ*, and leads to conjugate posterior computations
via the Gibbs sampler. The data-generating intensity is not assumed to be necessarily piecewise constant. Our methodology provides both a point estimate of the intensity function (posterior mean) and uncertainty quantification via marginal credible bands; see the figure below for an illustration.

<img src="https://raw.githubusercontent.com/mschauer/PointProcessInference.jl/master/assets/coal.png" alt="Intensity estimate for the mining data." width="67%">
*Illustration: Intensity estimation for the UK coal mining disasters data (1851-1962). The data are displayed via the rug plot in the upper margin of the plot, the posterior mean is given by a solid black line, while a 95% marginal credible band is shaded in light blue.*

## Installation

Change julia into the package manager mode by hitting `]`. Then run the command `add https://github.com/mschauer/PointProcessInference.jl`.

```
pkg> add https://github.com/mschauer/PointProcessInference.jl
```

## Usage

In the following example we load the UK coal mining disasters data and
performs its statistical analysis via the Poisson point process.


```
using PointProcessInference
using Random

Random.seed!(1234) # set RNG

observations, parameters, λinfo = PointProcessInference.loadexample("coal")

res = PointProcessInference.inference(observations; parameters...)
```

The main procedure has signature

```julia
ppinference(observations; title = "Poisson process", T = 1.0, n = 1, ...)
```

where `observations` is the sorted vector of Poisson event times, `T` is the endpoint of the time interval considered, and if
`observations` is an aggregate of `n` different independent observations (say aggregated for `n` days), this can be indicated by the parameter `n > 1`.

## High-quality plots

For high quality plotting, the package has a script `process-output-simple.jl` in the `contrib` folder that visualizes
the results with the help of `R` and `ggplot2`.
After installing the additional dependencies
```
pkg> add RCall
pkg> add DataFrames
```
call
```
include(joinpath(dirname(pathof(PointProcessInference)), "..", "contrib", "process-output-simple.jl"))
```
The script reads the variable `res`.


