
################ Generate or read data
function inference(observations_;
    title = "Poisson process", # optional caption for mcmc run
    summaryfile = nothing, # path to summaryfile or nothing
    T0 = 0.0, # start time
    T = maximum(observations), # end time
    n = 1, # number of aggregated samples in `observations`
    N = min(length(observations)÷4, 50), # number of bins
    IT = 30000, # number of iterations
    α1 = 0.1, β1 = 0.1, # parameters for Gamma Markov chain
    Π = Exponential(10), # prior on alpha
    τ = 0.7, # Set scale for random walk update on log(α)
    αind = 0.1, βind = 0.1, # parameters for the independence prior
    emp_bayes = false, # estimate βind using empirical Bayes
    verbose = true
)


    if T0 == 0
        observations = observations_
    else
        T -= T0
        observations = observations_ - T0
    end



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
    tt = @elapsed for i in 1:(IT-1)
        αψ = αζ = α[i]
        ζ[i,:] = updateζ(ψ[i,:], αψ, αζ)
        ψ[i+1,:] = updateψ(H, Δ, n, ζ[i,:], αψ, αζ, α1, β1)
        α[i+1], acc[i] = updateα(α[i], ψ[i + 1,:], ζ[i,:], τ,  Π)
    end

    if verbose
        println("Running time: $tt")
        println("")
        println("Average acceptance probability for updating: ",
            round(mean(acc); digits=3),"\n")
    end

    if summaryfile != nothing
        facc = open(summaryfile, "w")
        write(facc, "Data: ", string(title), "\n")
        write(facc, "Average acceptance probability: ", string(round(mean(acc); digits=3)), "\n")
        write(facc, "[T0, T, n, N] = ", string([T0, T, n, N]), "\n")
        write(facc, "Total number of events: ", string(sum(H)), "\n")
        write(facc, "tau = ", string(τ),"\n\n")
        write(facc, "Prior specification:", "\n")
        write(facc, "\talpha_ind = ", string(αind), "\n")
        write(facc, "\tbeta_ind = ", string(βind, "\n"))
        write(facc, "\talpha1 = ", string(α1), "\n")
        write(facc, "\tbeta1 = ", string(β1), "\n")
        write(facc, "\tPi = ", string(Π), "\n")
        close(facc)
    end

    return (title=title, observations = observations_, ψ = ψ, N = N, T0=T0, T = T + T0, breaks = breaks + T0, acc = acc)
end
