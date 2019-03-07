using DataFrames
using Statistics
using RCall

@rlibrary ggplot2

function plotposterior(res; figtitle="Markov chain prior", offset=0, hortics=0)
        ψ= res.ψ
        p = 0.05
        A = view(ψ, size(ψ, 1)÷2:size(ψ, 1), :)
        upper = mapslices(v -> quantile(v, 1 - p/2), A, dims=1)
        med = median(A, dims=1)
        ave = mean(A, dims=1)
        lower = mapslices(v -> quantile(v,  p/2), A, dims=1)
        breaks = offset.+ res.breaks
        N = length(breaks) - 1

        summaryψ=[ave' lower' med' upper']
        # dMarkov sets the bands of the posterior
        dMarkov = DataFrame(average=vec(ave), lower=vec(lower),
                            median=vec(med),  upper=vec(upper),
                            xmin=breaks[1:N], xmax=breaks[2:end])

        # dTrue contains the true intensity, evaluated on a grid
#        dTrue = DataFrame(x=gr, intensity=vals)

        # tMarkov sets the posterior mean (need to duplicate the final value)
        tMarkov = DataFrame(x=breaks, average=[vec(ave); vec(ave)[end]])
        obs = DataFrame(x=res.observations)

        #@rlibrary ggplot2  # cannot be called within the function appararently
        p = ggplot() + geom_rect(data=dMarkov,
        aes(xmin=:xmin, xmax=:xmax, ymin = :lower, ymax = :upper), fill = "lightsteelblue1") +
        geom_step(data=tMarkov, aes(x=:x, y=:average), colour="black", size=1.3) +
        ggtitle(figtitle) +
        xlab("") + ylab("") +
        geom_rug(data=obs, mapping=aes(x=:x), color="black", sides="t")
        if hortics == 0
            pp = p+ scale_x_continuous(limits=[minimum(breaks), maximum(breaks)])
        else
            pp = p + scale_x_continuous(limits=[minimum(breaks), maximum(breaks)], breaks=hortics)
        end
        pp
end
