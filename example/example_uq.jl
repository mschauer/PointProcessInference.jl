using Distributions
using PointProcessInference
using Plots
const PPI = PointProcessInference
include(PPI.plotscript()) #  basic plotting function

observations, parameters, λinfo = PointProcessInference.loadexample("testdat_n4000")
T = parameters.T
n = parameters.n

N = 1500
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N, τ=0.1)
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N, τ=0.1,Π=Exponential(100000))
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N, τ=0.1,Π=InverseGamma(.1))

plotposterior(res,λ=λinfo.λ)

# verify converges of α-chain
pyplot()
plot(res.α, label="α")


if false # testing with hdi intervals

    # test for HPD intervals
    using RCall
    @rlibrary HDInterval

    offset = 0
    ψ = res.ψ
    p = 0.05
    A = view(ψ, size(ψ, 1)÷2:size(ψ, 1), :)
    upper = mapslices(v-> quantile(v, 1 - p/2), A, dims=1)
    med = median(A, dims=1)
    ave = mean(A,dims=1)
    lower = mapslices(v-> quantile(v,  p/2), A, dims=1)
    breaks = offset.+ res.breaks
    N = length(breaks)-1

    @rput p
    @rput A
    R"""
    print(class(A))
    library(HDInterval)
    out = apply(A,2, function(x) {hdi(x, credMass=1-p)})
    """
    @rget out
    lower = out[1,:]'
    upper = out[2,:]'


    summaryψ=[ave' lower' med' upper']
    # dMarkov sets the bands of the posterior
    dMarkov = DataFrame(average=vec(ave), lower=vec(lower),
                        median=vec(med), upper=vec(upper),
                        xmin=breaks[1:N],xmax=breaks[2:end])
     λ=λinfo.λ
end
