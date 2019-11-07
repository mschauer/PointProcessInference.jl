using DelimitedFiles

workdir = @__DIR__
println(workdir)
cd(workdir)

# Make small test example
T = 10.0   # observation interval
n = 1   # number of copies of PPP observed
# Specify intensity function
λ = x ->  2* (5 + 4*cos(x)) * exp(-x/5)
λmax = λ(0.0)
observations =  PointProcessInference.samplepoisson((x,n)->λ(x)*n, n*λmax, 0.0, T, n)
sorted = false

n = 1
T = 10.0
observations = vec(readdlm("testdat_n1.csv"))

# At each iteration we keep track of a State
struct State
    modelindex::Int64
    logtarget::Float64
    ψ::Vector{Float64}
end

# precompute Δ en H for all models considered (could also do this 'on the fly', but that would amount to recomputing the same quantities many times)
Nmax = 40
Δvec = Vector{Float64}[]
Hvec = Vector{Int64}[]
for N in 1:Nmax
	breaks = range(0,T,length=N+1)
	push!(Δvec,diff(breaks))
	if sorted==true
	  push!(Hvec, PointProcessInference.counts_sorted(observations, breaks))
	else
	  push!(Hvec, PointProcessInference.counts(observations, breaks))
	end
end

"""
	Sample ψ from its posterior distribution within a particular model

	(α,β): are parameters of Gamma-prior on bin size
	n: nr  of observations
	H, Δ:  bin characteristics for posterior pars
"""
postψ(H,Δ,α,β,n) = [rand(Gamma(α+H[k], 1.0/(β+n*Δ[k]))) for k in eachindex(H)]


"""
	log π(N), presently all models are equally likely
"""
function logpriormodel(N)
	0.0
end

"""
	proposalratio when proposing move from N to Nᵒ

	Computes
```math
\\frac{q(N \\mid N^\\circ)}{q(N^\\circ \\mid N)}
```
"""
function proposalratio(N, Nᵒ;η=0.4)
	if N==1 & Nᵒ==2
		return(2η)
	elseif N==2 & Nᵒ==1
		return(0.5/η)
	else
		return(1.0)
	end
end

"""
	At N=1, stay with prob 0.5 in 1, else move to 2.
	At all models N>=2, with probability η move to N+1, with probability η move to N-1, with probability 1-2η stay at N
"""
function modelindexproposal(N;η=0.4)
	u = rand()
	if N==1
		return(ifelse(u<0.5,1,2))
	else
		if u < η
			return(N-1)
		elseif u > 1-η
			return(N+1)
		else
			return(N)
		end
	end
end

α = β = 0.1
Ninit = 2
logtargetinit = PointProcessInference.mloglikelihood(Ninit, observations,T, n, α, β) + logpriormodel(Ninit)
ψinit = postψ(Hvec[Ninit],Δvec[Ninit],α,β,n)
states = [State(Ninit,logtargetinit,ψinit)]

ITER = 1000
η = 0.45 # prob of moving to other model

breaksvec = Float64[]
ψvec = Float64[]
itervec = Int64[]

for i in 2:ITER
	global breaksvec
	global ψvec
	global itervec
	N = states[i-1].modelindex
	Nᵒ = modelindexproposal(N;η=η)
	print("propose ", N, " to ", Nᵒ)
	logtargetᵒ = PointProcessInference.mloglikelihood(Nᵒ, observations,T, n, α, β) + logpriormodel(Nᵒ)
	A = logtargetᵒ - states[i-1].logtarget + log(proposalratio(N,Nᵒ;η=η))
	if (log(rand())<A)
		ψᵒ = postψ(Hvec[Nᵒ],Δvec[Nᵒ],α,β,n)
		push!(states, State(Nᵒ,logtargetᵒ,ψᵒ))
		println("   acc")
	else
		ψ = postψ(Hvec[N],Δvec[N],α,β,n)
		push!(states, State(N,states[i-1].logtarget,ψ))
		println("   --")
	end

	St = states[i]
	breaksvec = vcat(breaksvec, collect(range(0,T,length=St.modelindex+1)))
	ψvec = vcat(ψvec, vcat(St.ψ,St.ψ[end]))
	itervec = vcat(itervec, fill(i,St.modelindex+1))


end




 # plot a particular State
using DataFrames
using RCall

# St = states[35]
# df = DataFrame(x=collect(range(0,T,length=St.modelindex+1)),y=vcat(St.ψ,St.ψ[end]))

df = DataFrame(x=breaksvec, y= ψvec, iter=itervec)
@rput df



# extract states during iterations
dfmodel = DataFrame(model= map(x->x.modelindex, states), iter=collect(1:ITER))
@rput dfmodel

# gr = collect(0:0.01:10)
# dftrueintensity = DataFrame(x=gr,y=λ.(gr) )
# @rput dftrueintensity

mloglik = [PointProcessInference.mloglikelihood(Nᵒ, observations,T, n, α, β) for Nᵒ in 1:Nmax]
mll_df = DataFrame(x=1:Nmax,y=mloglik)
@rput mll_df

R"""
library(tidyverse)
library(ggplot2)
library(gridExtra)

trueintens <- function(x) {
	2* (5 + 4*cos(x)) * exp(-x/5)
}


p1 <- df  %>%
mutate(iteration=as.factor(iter)) %>%
ggplot(aes(x=x,y=y,group=iteration,colour=iter)) + geom_step() +
 stat_function(fun = trueintens,colour='red')

p2 <- df %>% dplyr::filter(iter>500) %>%
mutate(iteration=as.factor(iter)) %>%
ggplot(aes(x=x,y=y,group=iteration,colour=iter)) + geom_step() +
 stat_function(fun = trueintens,colour='red')

p3 <- dfmodel %>% ggplot(aes(x=iter,y=model)) + geom_point()

p4 <- mll_df %>% ggplot(aes(x=x,y=y)) + geom_point() + ggtitle("marginal loglikelihood")

grid.arrange(p1,p2,p3,p4,ncol=1)
"""
