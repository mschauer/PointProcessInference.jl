"""
   samplepoisson

Sample a non homogeneous Poisson process on [0,T] via thinning.
λ is the intensity;
λmax is an upper bound for λ(x), when x in [0,T].
args contains additional arguments passed to the function λ
"""
function samplepoisson(λ::Function, λmax, T, args...)
    t = 0.0
    tt = zeros(0)
    while t <= T
        t = t - log(rand())/λmax
        if rand() ≤ λ(t, args...)/λmax
            push!(tt, t)
        end
    end
    tt
end



"""
    counts(xx, grid)

Count how many points fall between the grid points in `grid`.

# Example:

```
  julia> counts(rand(10), 0:0.5:1)
  2-element Array{Int64,1}:
 6
 4
```
"""
function counts(xx, grid)
    c = zeros(Int, length(grid) + 1)
    for x in xx
        c[first(searchsorted(grid, x))] += 1
    end
    c[2:length(grid)]
end


function counts_sorted(xx, grid)
    n = length(grid)
    c = zeros(Int, n + 1)
    i = 1
    for x in xx
        while i <= n && x > grid[i]
            i += 1
        end
        c[i] += 1
    end
    c[2:length(grid)]
 end

################# functions for updating

"""
    updateψ

Sample from the distribution of ψ, conditional on ζ, αψ and αζ.

Arguments:
- H::Array{Int64}      (count of events over bins)
- Δ::Array{Float64}    (lengths of bins)
- n::Int64             (sample size)
- ζ::Array{Float64}
- αψ::Float64
- αζ::Float64
- α1::Float64          (shape parameter of Gamma prior on ψ[1])
- β1::Float64          (rate parameter of Gamma prior on ψ[1])
"""
function updateψ(H::Array{Int64}, Δ::Array{Float64}, n::Int64, ζ::Array{Float64},
    αψ::Float64, αζ::Float64, α1::Float64, β1::Float64
)
    N = length(H)
    ψ = zeros(N)
    ψ[1] = rand(Gamma(α1 + αζ + H[1], 1.0/(β1 + αζ/ζ[2] + n*Δ[1])))
    for k in 2:(N-1)
        ψ[k] = rand(Gamma(αψ + αζ + H[k], 1.0/(αψ/ζ[k] + αζ/ζ[k+1] + n*Δ[k])))
    end
    ψ[N] = rand(Gamma(αψ + H[N], 1.0/(αψ/ζ[N] + n*Δ[N])))
    ψ
end

"""
    updateζ

Sample from the distribution of ζ, conditional on ψ, αψ and αζ.

Arguments:
- ψ::Array{Float64}
- αψ::Float64
- αζ::Float64
"""
function updateζ(ψ::Array{Float64}, αψ::Float64, αζ::Float64)
    N = length(ψ)
    ζ = zeros(N)
    for k in 2:N
        ζ[k] = rand(InverseGamma(αψ + αζ, αζ*ψ[k-1] + αψ*ψ[k]))
    end
    ζ
end


"""
    sumψζ

Helper function for updating α

Computes
```math
\\sum_{k=2}^N \\log \\frac{ψ_{k-1}ψ_k)}{ζ_k^2} -\\frac{ψ_{k-1}ψ_k}{ζ_k}
```

Arguments:
- ψ::Array{Float64}
- ζ::Array{Float64}
"""
function sumψζ(ψ::Array{Float64},ζ::Array{Float64})
    res = 0.0
    for k in 2:length(ψ)
        res += log(ψ[k-1]) + log(ψ[k]) - 2*log(ζ[k]) - (ψ[k-1] + ψ[k])/ζ[k]
    end
    res
end

function log_q(α::Float64, N::Int64, sumval::Float64)
    2*(N - 1)*(α*log(α) - lgamma(α)) + α*sumval
end

function log_qtilde(lα::Float64, N::Int64, sumval::Float64, Π)
    log_q(exp(lα), N, sumval) + lα + logpdf(Π, exp(lα))
end

# given ψ and ζ, draw α using symmetric random walk on log(α)
function updateα(α::Float64, ψ::Array{Float64}, ζ::Array{Float64}, τ::Float64, Π)
    sumval = sumψζ(ψ, ζ)
    N = length(ψ)
    lα = log(α)
    ll = log_qtilde(lα, N, sumval, Π)
    lα_prop = lα + τ * randn()
    ll_prop = log_qtilde(lα_prop, N, sumval, Π)
    if log(rand()) < (ll_prop - ll)
        return exp(lα_prop), true
    else
        return α, false
    end
end

# compute marginal likelihood for N=2..Nmax, here Nmax should be >=2
function marginal_loglikelihood(Nmax::Int64, observations::Vector{Float64},
                       T::Float64, n::Int64, αind::Float64, βind::Float64)
    mll = zeros(Nmax-1)
    for N in 2:Nmax
        breaks = range(0, stop=T, length=N+1)
        Δ = diff(breaks)
        H = counts(observations, breaks)
        ltip = lgamma.(αind .+ H) .- (αind .+ H) .* log.(n*Δ .+ βind) # ltip = log terms in product
        mll[N-1] = T*n + αind*N*log(βind) - N*lgamma(αind) + sum(ltip)
    end
    2:Nmax, mll
end

function elpd_DIC(Nmax::Int64, observations::Vector{Float64},
        T::Float64, n::Int64, αind::Float64, βind::Float64)
    elpd = Vector{Float64}(Nmax-1)
    for N in 2:Nmax
        breaks = range(0, stop=T, length=N+1)
        Δ = diff(breaks)
        H = counts(observations, breaks)
        tip = (αind + H)./(n*Δ + βind)  # tip = terms in product
        ll = sum(H.*log.(tip) - n*Δ .* tip)
        νDIC = 2 * sum(H.*(log.(αind + H) - digamma.(αind + H)) )
        elpd[N-1] = ll - νDIC
    end
    2:Nmax, elpd
end

# Determine N as the largest N for which each bin has at least Nmin observations
function determine_number_breaks(Nmin::Int64, T::Float64, observations::Vector{Float64})
    too_many = true
    N = 1
    while too_many
        breaks = range(0, stop=T, length=N+1) # determine breaks points with N bins
        H = counts(observations, breaks)  # compute counts over bins
        if minimum(H) >= Nmin
            N += 1 # we may be able to use more bins
        else
            too_many = false # each bin contains at least Nmin observations
            return N - 1
        end
    end
end

"""
    ebβ

Determine β by maximising the marginal likelihood, for fixed α and N.
"""
function ebβ(α,  H, Δ, n, N)
    GG(α, H, Δ, n, N) = (lβ) -> -α*N*lβ[1] + N * lgamma(α) -  sum(lgamma.(α+H)) + sum((α+H).* log.(n*Δ+exp(lβ[1])))
    result = optimize(GG(α, H, Δ, n, N), [0.0], BFGS())
    exp.(result.minimizer)[1]
end

"""
    ebα

Determine α by maximising the marginal likelihood, for fixed β and N.
"""
function ebα(β,  H, Δ, n, N)
    GG(β, H, Δ, n, N) = (lα) -> -exp(lα[1])*N* log(β) + N * lgamma(exp(lα[1])) -  sum(lgamma.(exp(lα[1])+H)) + sum((exp(lα[1])+H).* log.(n*Δ+β))
    result = optimize(GG(β, H, Δ, n, N), [0.0], BFGS())
    exp.(result.minimizer)[1]
end
