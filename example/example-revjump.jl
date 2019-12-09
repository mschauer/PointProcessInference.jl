 #add https://github.com/mschauer/PointProcessInference.jl#revjump
using DelimitedFiles
using PointProcessInference
const PPI = PointProcessInference
using Distributions
using DataFrames
using RCall

workdir = @__DIR__
println(workdir)
cd(workdir)

observations, parameters, λinfo = PointProcessInference.loadexample("testdat_n1")  # "testdat_n1"
T = parameters.T
n = parameters.n

########## run the rev jump algorithm
αind = 0.1; βind = 0.1
Nmax = 50
priorN = Poisson(23)#DiscreteUniform(1,Nmax) #
ITER = 30_000
T0 = 0.0
states, df = revjump(observations,T0, T,n,priorN; ITER=ITER, αind=αind, βind=βind,Nmax=Nmax)


########## make dataframes for plotting
# df for models visited
dfmodel = DataFrame(model= map(x->x.modelindex, states), iter=collect(1:ITER))

# df for marginal likelihood
mloglik = [PPI.mloglikelihood(N, observations,T0, T, n, αind, βind) for N in 1:Nmax]
mll_df = DataFrame(x=1:Nmax,y=mloglik)

# df for posterior mean and quantiles
Nfinest = maximum(dfmodel.model)  # finest level
xgrid = collect(range(0,T,step=0.001))
Ngrid = length(xgrid)
BURNIN = div(ITER,2)
out = zeros(ITER-BURNIN,Ngrid)  # rows infex iterations
for i in BURNIN+1:ITER
	out[i-BURNIN,:] = evalstepfunction(states[i].ψ,T).(xgrid)
end
p = 0.05
upper = vec(mapslices(v-> quantile(v, 1 - p/2), out; dims= 1))
ave = vec(mean(out,dims=1))
lower = vec(mapslices(v-> quantile(v,  p/2), out; dims= 1))
dout = DataFrame(lower=lower,upper=upper, ave=ave,x=xgrid)

# df for observations
obsdf = DataFrame(x=observations)

########## R based plotting:
@rput dfmodel
@rput df
@rput mll_df
@rput dout
@rput T
@rput obsdf
@rput dfmodel
@rput Nfinest

R"""
library(tidyverse)
library(ggplot2)
library(gridExtra)
theme_set(theme_minimal())

trueintens <- function(x) {
	2* (5 + 4*cos(x)) * exp(-x/5)
}

# plot of iterations
#	p1 <- df  %>%
#	mutate(iteration=as.factor(iter)) %>%
#	ggplot(aes(x=x,y=y,group=iteration,colour=iter)) + geom_step() +
#	 stat_function(fun = trueintens,colour='red')

#	p2 <- df %>% dplyr::filter(iter>500) %>%
#	mutate(iteration=as.factor(iter)) %>%
#	ggplot(aes(x=x,y=y,group=iteration,colour=iter)) + geom_step() +
#	 stat_function(fun = trueintens,colour='red') */

p4 <- mll_df %>% ggplot(aes(x=x,y=y)) + geom_path() +
  scale_x_continuous(name="N")+ ylab("marginal loglikelihood")

pdf('marginal_loglik.pdf',width=8,height=4)
	show(p4)
dev.off()

# plot posterior mean and quantiles
p <- dout %>%  ggplot() +
	geom_ribbon(aes(x=x,ymin=lower,ymax=upper),fill="lightsteelblue1") +
	geom_path(aes(x=x,y=ave),colour='black',size=1.3) +
	 stat_function(fun = trueintens,colour='red',size=1.5)+
		#ggtitle("Discrete Uniform{1,...,50}")+
		ggtitle("ShiftPoisson(23)")+
	   theme(plot.title = element_text(hjust = 0.5))+xlab("")+ylab("")+
xlim(0,T-10^(-5)) +  geom_rug(data=obsdf, mapping=aes(x=x), color="black",sides='t')

pdf('estimates_revjump.pdf',width=4,height=4)
	show(p)
dev.off()

# plot models visited over iterations
p5 <- dfmodel %>% ggplot(aes(x=iter,y=model)) + geom_path() +
 xlab("Iteration") + scale_y_continuous(name="N")

pdf('estimates_revjump_modelsvisited.pdf',width=8,height=4)
 	show(p5)
dev.off()

# plot relative freq of models visited over iterations
p6 <- dfmodel %>% mutate(modelf=as.factor(model)) %>% ggplot(aes(x=modelf))+
   geom_bar(aes(y = (..count..)/sum(..count..))) +xlab("N") + ylab("Rel. frequency")

pdf('estimates_revjump_modelsvisited_bar.pdf',width=8,height=2)
	show(p6)
dev.off()
"""
