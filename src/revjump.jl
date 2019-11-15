# At each iteration we keep track of a State
struct State
    modelindex::Int64
    logtarget::Float64
    ψ::Vector{Float64}
end

"""
 	precompute Δ en H for all models considered (could also do this 'on the fly', but that would amount to recomputing the same quantities many times)
"""
function computebinning(T; Nmax = 40, sorted=false)
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
	Δvec, Hvec
end


"""
	Sample ψ from its posterior distribution within a particular model

	(α,β): are parameters of Gamma-prior on bin size
	n: nr  of observations
	H, Δ:  bin characteristics for posterior pars
"""
postψ(H,Δ,α,β,n) = [rand(Gamma(α+H[k], 1.0/(β+n*Δ[k]))) for k in eachindex(H)]


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
function modelindexproposal(N; η=0.4)
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

"""
	observations: vector of observed times
"""
function revjump(observations,T,n, Hvec, Δvec, priorN; ITER=10000, Ninit =2, α = 0.1, β = 0.1, η=0.45)
	logtargetinit = PointProcessInference.mloglikelihood(Ninit, observations,T, n, α, β) +
							logpdf(priorN,Ninit)
	ψinit = postψ(Hvec[Ninit],Δvec[Ninit],α,β,n)
	states = [State(Ninit,logtargetinit,ψinit)]

	breaksvec = Float64[]
	ψvec = Float64[]
	itervec = Int64[]

	for i in 2:ITER
		N = states[i-1].modelindex
		Nᵒ = modelindexproposal(N; η=η)
		logtargetᵒ = PointProcessInference.mloglikelihood(Nᵒ, observations,T, n, α, β) +
							logpdf(priorN,Nᵒ)
		A = logtargetᵒ - states[i-1].logtarget + log(proposalratio(N,Nᵒ;η=η))
		if (log(rand())<A)
			ψᵒ = postψ(Hvec[Nᵒ],Δvec[Nᵒ],α,β,n)
			push!(states, State(Nᵒ,logtargetᵒ,ψᵒ))
		else
			ψ = postψ(Hvec[N],Δvec[N],α,β,n)
			push!(states, State(N,states[i-1].logtarget,ψ))
		end
		St = states[i]
		breaksvec = vcat(breaksvec, collect(range(0,T,length=St.modelindex+1)))
		ψvec = vcat(ψvec, vcat(St.ψ,St.ψ[end]))
		itervec = vcat(itervec, fill(i,St.modelindex+1))
	end
	states, DataFrame(x=breaksvec, y= ψvec, iter=itervec)
end

# evaluate stepfunction with weights ψ at point x ∈ [0,T], assuming length(ψ) equally sized bins
function evalstepfunction(x,ψ,T)
	if (x<0.0) | (x>=T)
		return(0.0)
	else
		widthbin = T/length(ψ)
		binindex = Int(div(x,widthbin) +1)
		return(ψ[binindex])
	end
end
evalstepfunction(ψ,T) = (x) -> evalstepfunction(x,ψ,T)
