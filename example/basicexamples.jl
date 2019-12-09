using Distributions
using PointProcessInference
const PPI = PointProcessInference

# example with generated data
obs = rand(Exponential(2), 100)
res = PPI.inference(obs)

include(PPI.plotscript())
plotposterior(res)

# example from the accompanying paper, to see available examples type "?PPI.loadexample"

observations, parameters, λinfo = PPI.loadexample("testdat_n1")
N = 15 # nr of bins
res = PPI.inference(observations, title = "Poisson process", T = parameters.T, n = parameters.n, N=N, samples = 1:1:10000, τ=0.6)

include(PPI.plotscript())
plotposterior(res,λ=λinfo.λ)
