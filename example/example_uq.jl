using Distributions
using PointProcessInference
const PPI = PointProcessInference

observations, parameters, λinfo = PointProcessInference.loadexample("testdat_n4000")
T = parameters.T
n = parameters.n

res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=200)

include(PointProcessInference.plotscript())
plotposterior(res)
