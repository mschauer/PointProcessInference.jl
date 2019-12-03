using Distributions
using PointProcessInference
using Plots
const PPI = PointProcessInference
using DataFrames
using Statistics
using RCall
@rlibrary ggplot2

workdir = @__DIR__
println(workdir)
cd(workdir)


# make specific function for plotting the results in this experiment
function plotposterior2(res;  figtitle="Markov chain prior", λ=0, p=0.05)
        ψ = res.ψ
        A = view(ψ, size(ψ, 1)÷2:size(ψ, 1), :)
        upper = mapslices(v-> quantile(v, 1 - p/2), A, dims=1)
        med = median(A, dims=1)
        ave = mean(A,dims=1)
        lower = mapslices(v-> quantile(v,  p/2), A, dims=1)
        breaks = offset.+ res.breaks
        N = length(breaks)-1

        summaryψ=[ave' lower' med' upper']
        # dMarkov sets the bands of the posterior
        dMarkov = DataFrame(average=vec(ave), lower=vec(lower),
                            median=vec(med), upper=vec(upper),
                            xmin=breaks[1:N],xmax=breaks[2:end])
        # dTrue contains the true intensity, evaluated on a grid

        if !(λ==0)
            gr = collect(range(minimum(breaks),step=0.02,stop=maximum(breaks)))
            dTrue = DataFrame(x=gr, intensity=λ.(gr .- offset))
        end
        # tMarkov sets the posterior mean (need to duplicate the final value)
        tMarkov = DataFrame(x=breaks, average=[vec(ave); vec(ave)[end]])
        obs = DataFrame(x=res.observations .+ offset)
        # make basic plot
        @rput dMarkov
        @rput tMarkov
        @rput dTrue
        @rput figtitle
        R"""
        library(ggplot2)
        theme_set(theme_minimal())
        p = ggplot() + geom_rect(data=dMarkov,
          aes(xmin=xmin,xmax=xmax,ymin = lower, ymax = upper), fill = "lightsteelblue1") +
          geom_step(data=tMarkov, aes(x=x,y=average),colour="black",size=0.3)+
           theme(plot.title = element_text(hjust = 0.5)) +
          xlab("")+ylab("")+ geom_line(data=dTrue, aes(x=x, y=intensity),size=.4,colour="red")+
            ggtitle(figtitle)
          #geom_rug(data=obs, mapping=aes(x=:x), color="black",sides="t")+theme_light()+


        pdf("estimates.pdf",width=4,height=4)
             show(p)
         dev.off()
         """

end



observations, parameters, λinfo = PointProcessInference.loadexample("testdat_n4000")
T = parameters.T
n = parameters.n

N = 500
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N, τ=0.1)
plotposterior2(res,λ=λinfo.λ, figtitle="N=500, Exp(0.1)")

N = 1000
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N, τ=0.1)
plotposterior2(res,λ=λinfo.λ, figtitle="N=1000, Expl(0.1)")

N = 500
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N, τ=0.1,Π=InverseGamma(.1))
plotposterior2(res,λ=λinfo.λ, figtitle="N=500, IG(0.1,1)")

N = 1000
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N, τ=0.1,Π=InverseGamma(.1))
plotposterior2(res,λ=λinfo.λ, figtitle="N=1000, IG(0.1,1)")

N = 500
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N, τ=0.1,Π=Levy())
plotposterior2(res,λ=λinfo.λ, figtitle="N=500, Levy(0,1)")

N = 1000
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N, τ=0.1,Π=Levy())
plotposterior2(res,λ=λinfo.λ, figtitle="N=1000, Levy(0,1)")

N = 5000
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N, τ=0.05,Π=Levy())
plotposterior2(res,λ=λinfo.λ, figtitle="N=5000, Levy(0,1)")

################
observations, parameters, λinfo = PointProcessInference.loadexample("testdat_n1")
T = parameters.T
n = parameters.n
N=30
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N, τ=0.5,Π=Levy(0,1))
plotposterior2(res,λ=λinfo.λ, figtitle="N=30, Levy(0,1)")

N=100
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N, τ=0.5,Π=Levy(0,1))
plotposterior2(res,λ=λinfo.λ, figtitle="N=100, Levy(0,1)")

N=30
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N)
plotposterior2(res,λ=λinfo.λ, figtitle="N=30, Exp(0.1)")

N=100
res = PPI.inference(observations, title = "Poisson process", T = T, n = n, N=N)
plotposterior2(res,λ=λinfo.λ, figtitle="N=100, Exp(0.1)")


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
