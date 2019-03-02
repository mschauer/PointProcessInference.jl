to = TimerOutput()

Random.seed!(1234) # set RNG

################ Generate or read data
function ppp(observations;
    title = "Poisson process", # optional caption for mcmc run
    summaryfile = nothing, # path to summaryfile or nothing
    T = maximum(observations), # total time
    n = 1, # number of aggregated samples in `observations`
    N = min(length(observations)÷4, 50), # number of bins
    IT = 30000, # number of iterations
    α1 = 0.1, β1 = 0.1, # parameters for Gamma Markov chain
    Π = Exponential(10), # prior on alpha
    τ = 0.7, # Set scale for random walk update on log(α)
    αind = 0.1, βind = 0.1, # parameters for the independence prior
    emp_bayes = false # estimate βind using empirical Bayes
)


    ################ Data processing

    breaks = range(0, stop=T, length=N+1) # linspace(0,T,N+1)
    Δ = diff(breaks)

    # if the observations are sorted, the bin counts can be computed faster
    if issorted(observations)
        sorted = true
        H = counts_sorted(observations, breaks)    # extract vector H
    else
        sorted = false
        H = counts(observations, breaks)           # extract vector H
    end

    if emp_bayes == true
        βind = ebβ(αind, H, Δ, n, N)
    end


    ################## Specification number of bins N

    # option 1a: maximise marginal log-likelihood with independence prior
    Nmax = length(observations)÷2
    Nvals, mll = marginal_loglikelihood(Nmax, observations, T, n, αind, βind)
    Nopt = Nvals[argmax(mll)]


    ################### Initialisation of algorithms

    # nr of iterations
    ψ = zeros(IT, N)  # each row is an iteration
    ζ = zeros(IT - 1, N)
    α = zeros(IT)
    α[1] = 1.0  # initial value for α
    acc = zeros(Bool, IT - 1)  # keep track of MH acceptance

    # Initialise, by drawing under independence prior
    post_ind = zeros(N, 2)
    for k in 1:N
    	post_ind[k,:] = [αind + H[k], βind + n*Δ[k]]  # note that second parameter is the rate
    	ψ[1,k] = rand(Gamma(post_ind[k,1], 1.0/(post_ind[k,2])))
    end

    # Gibbs sampler
    for i in 1:(IT-1)
        αψ = αζ = α[i]
        @timeit to "update ζ"  ζ[i,:] = updateζ(ψ[i,:], αψ, αζ)
        @timeit to "update ψ"  ψ[i+1,:] = updateψ(H, Δ, n, ζ[i,:], αψ, αζ, α1, β1)
        @timeit to "update α"  α[i+1], acc[i] = updateα(α[i], ψ[i + 1,:], ζ[i,:], τ,  Π)
    end

    println("Timing mcmc steps:")
    show(to, allocations = false, compact = true)

    println("")
    println("Average acceptance probability for updating  equals: ",
    round(mean(acc);digits=3),"\n")

    if summaryfile != nothing
        facc = open(summaryfile,"w")
        write(facc, "data: ", string(title),"\n")
        write(facc, "Average acceptance probability equals: ",string(round(mean(acc);digits=3)),"\n")
        write(facc, "[T, n,  N] = ",string([T, n, N]),"\n")
        write(facc, "total number of events ",string(sum(H)),"\n")
        write(facc, "tau = ",string(τ),"\n\n")
        #write(facc, "elapsed time ",string(elapsed_time), "\n\n")
        write(facc, "---- Prior specification ----","\n")
        write(facc, "alpha_ind = ",string(αind),"\n")
        write(facc, "beta_ind = ",string(βind,"\n"))
        write(facc, "alpha1 = ",string(α1),"\n")
        write(facc, "beta1 = ",string(β1),"\n")
        write(facc, "Pi = ",string(Π),"\n")
        close(facc)
    end

    return (ψ = ψ , acc = acc)
end
