using Distributions
using PointProcessInference
const PPI = PointProcessInference

observations, parameters, λinfo = PointProcessInference.loadexample("testdat_n4000")
T = parameters.T
n = parameters.n

N = 500
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N, τ=0.1)

include(PPI.plotscript())
plotposterior(res,λ=λinfo.λ)
