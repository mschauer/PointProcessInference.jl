to = TimerOutput()

Random.seed!(1234) # set RNG

################ Generate or read data
function ppp(observations; T = 1.0, n = 1, λ = nothing, λmax = NaN)

    emp_bayes = false
    sorted = false  # if the observations are sorted, the bin counts can be computed faster by setting this to true

    ################ Prior specification
    N = min(Int64(round(length(observations)/4)), 50)
    breaks = range(0,stop=T, length=N+1) # linspace(0,T,N+1)
    Δ = diff(breaks)
    if sorted==true
      H = counts_sorted(observations, breaks)           # extract vector H
    else
      H = counts(observations, breaks)           # extract vector H
    end

    αind = 0.1
    if emp_bayes == true
        βind = ebβ(αind,  H, Δ, n, N)
    else
        βind = 0.1
    end

    # Set parameters for Gamma Markov chain
    α1 = 0.1
    β1 = 0.1

    # prior on alpha
    Π = Exponential(10) # julia uses scale as parameter instead of rate

    ################# Set scale for random walk update on log(α)
    τ = 0.7
    ################## Specification number of bins N

    # option 1a: maximise marginal log-likelihood with independence prior
    Nmax = Int64(round(length(observations)/2;digits=0))
    Nvals, mll = marginal_loglikelihood(Nmax,observations, T, n, αind, βind)
    Nopt = Nvals[argmax(mll)]
    #println(Nopt)

    ################### Initialisation of algorithms


    IT = 30000 # nr of iterations
    ψ = zeros(IT,N)  # each row is an iteration
    ζ = zeros(IT-1,N)
    α = zeros(IT)
    α[1] = 1.0  # initial value for α
    acc = Vector{Int8}(undef,IT-1)  # keep track of MH acceptance

    # Initialise, by drawing under independence prior
    post_ind = zeros(N,2)
    for k=1:N
    	post_ind[k,:] = [αind+H[k], βind+n*Δ[k]]  # note that second parameter is the rate
    	ψ[1,k] = rand(Gamma(post_ind[k,1], 1.0/(post_ind[k,2])))
    end

    # Gibbs sampler
    for i=1:(IT-1)
      αψ = α[i];  αζ = α[i];
    @timeit to "update ζ"  ζ[i,:] = updateζ(ψ[i,:],αψ,αζ)
    @timeit to "update ψ"  ψ[i+1,:] = updateψ(H, Δ, n, ζ[i,:], αψ, αζ, α1, β1)
    @timeit to "update α"  α[i+1], acc[i] = updateα(α[i], ψ[i+1,:],ζ[i,:], τ,  Π)
    end

    println("Timing mcmc steps:")
    show(to, allocations = false, compact = true)

    println("")
    println("Average acceptance probability for updating  equals: ",
    round(mean(acc);digits=3),"\n")

    facc = open("./info.txt","w")
    write(facc, "data: ", string(data_choice),"\n")
    write(facc, "Average acceptance probability equals: ",string(round(mean(acc);digits=3)),"\n")
    write(facc, "[T, n,  N] = ",string([T, n, N]),"\n")
    write(facc, "total number of events ",string(sum(H)),"\n")
    write(facc, "tau = ",string(τ),"\n\n")
    #write(facc, "elapsed time ",string(elapsed_time), "\n\n")
    write(facc, "---- Prior specification ----","\n")
    write(facc, "alpha_ind= ",string(αind),"\n")
    write(facc, "beta_ind= ",string(βind,"\n"))
    write(facc, "alpha1= ",string(α1),"\n")
    write(facc, "beta1= ",string(β1),"\n")
    write(facc, "Pi= ",string(Π),"\n")
    close(facc)

    return (ψ = ψ , acc = acc)
end

# Process output
#include("process-output.jl")
include("process-output-simple.jl")
