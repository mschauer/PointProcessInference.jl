using DataFrames
using Statistics
using RCall

function showresults(res; p = .05)
    observations = res.observations
    breaks = res.breaks
    title = res.title
    ψ = res.ψ
    N = res.N
    acc = res.acc

    A = view(ψ, size(ψ, 1)÷2:size(ψ, 1), :)
    upper = mapslices(v-> quantile(v, 1 - p/2), A, dims=1)
    med = median(A, dims=1)
    ave = mean(A,dims=1)
    lower = mapslices(v-> quantile(v,  p/2), A, dims=1)

    ## plotting

    # dMarkov sets the bands of the posterior
    dMarkov = DataFrame(average=vec(ave), lower=vec(lower),
                        median=vec(med), upper=vec(upper),
                        xmin=breaks[1:N],xmax=breaks[2:end])
    # tMarkov sets the posterior mean (need to duplicate the final value)
    tMarkov = DataFrame(x=breaks, average=[vec(ave); vec(ave)[end]])

    obs = DataFrame(x=observations)

    @rlibrary ggplot2
    p = ggplot() + geom_rect(data=dMarkov,
    aes(xmin=:xmin,xmax=:xmax,ymin = :lower, ymax = :upper), fill = "lightsteelblue1") +
    geom_step(data=tMarkov, aes(x=:x,y=:average),colour="black",size=1.3)+
    ggtitle("Markov chain prior")+
    xlab("")+ylab("")+
    geom_rug(data=obs, mapping=aes(x=:x), color="black",sides="t")

    if data_choice in ["generated","testdat_n1","testdat_n5","testdat_n4000","simpson_n200","simpson_n500", "unit_intensity","unit_intensity_withpeak"]
        gr = range(0,stop=T,length=1000); vals = λ.(gr)
        dTrue = DataFrame(x=gr, intensity=vals)
        show(p + geom_line(data=dTrue, aes(x=:x,y=:intensity),colour="red",linetype = "solid",size=1.5) )
    else
        show(p)
    end
end
