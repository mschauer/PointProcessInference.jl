

workdir = @__DIR__
println(workdir)
cd(workdir)

T = 10.0   # observation interval
n = 50   # number of copies of PPP observed
# Specify intensity function
λ = x ->  2* (5 + 4*cos(x)) * exp(-x/5)
λmax = λ(0.0)
observations =  PointProcessInference.samplepoisson((x,n)->λ(x)*n, n*λmax, 0.0, T, n)
sorted = false

struct State
    modelindex::Int64
    logtarget::Float64
    ψ::Vector{Float64}
end

Nmax = 40

# precompute Δ en H for all models
Δvec = Vector{Float64}[]#zeros(Nmax)
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
	Sample ψ from its posterior distribution
"""
function postψ(H,Δ,α,β,n)
	N = length(H)
	out = zeros(N)
	for k in 1:N
		out[k] = rand(Gamma(α+H[k], 1.0/(β+n*Δ[k])))
	end
	out
end

"""
	log π(N)
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

ITER = 500
η = 0.45 # prob of moving to other model

for i in 2:ITER
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
		push!(states, states[i-1])
		println("   --")
	end
end
