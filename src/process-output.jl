print("Average acceptance probability equals: ",mean(acc),"\n")

############# Write output to csv files (ensure a relative directory with name "out" exists)

facc = open("./out/info.txt","w")
write(facc, "data: ", string(data_choice),"\n")
write(facc, "Average acceptance probability equals: ",string(round(mean(acc);digits=3)),"\n")
write(facc, "[T, n,  N] = ",string([T, n, N]),"\n")
write(facc, "total number of events ",string(sum(H)),"\n")
write(facc, "tau = ",string(τ),"\n\n")
write(facc, "elapsed time ",string(elapsed_time), "\n\n")
write(facc, "---- Prior specification ----","\n")
write(facc, "alpha_ind= ",string(αind),"\n")
write(facc, "beta_ind= ",string(βind,"\n"))
write(facc, "alpha1= ",string(α1),"\n")
write(facc, "beta1= ",string(β1),"\n")
write(facc, "Pi= ",string(Π),"\n")
close(facc)

if data_choice in ["generated","testdat_n1","testdat_n5","testdat_n4000","simpson_n200","simpson_n500", "unit_intensity","unit_intensity_withpeak"]
    gr = range(0,T,length=1000); vals = λ.(gr)
else
    gr = randn(3); vals = randn(3)
end

writedlm("./out/true_intensity.csv",[gr vals],',')   # write true intensity function on fine grid to csv
writedlm("./out/observations.csv",observations,',')   # observations
writedlm("./out/breaks.csv",breaks,',')      # breaks specifying the binning (including begin and end point)
writedlm("./out/post_ind.csv",post_ind,',')  # posterior parameter under independence prior
writedlm("./out/alpha.csv",α,',')            # all iterates from alpha
writedlm("./out/marglikelihood_N.csv",[Nvals mll],',')            # marginal likelihood values


### Randomly sample 3 values of psi and zeta and write these to csv files
indices = sort(sample(2:N,3,replace=false))

fψ = open("./out/psi.csv","w")
head = prod("psi$i"*(i == indices[end] ? "\n" : ", ")  for i in indices)
write(fψ, head)
writedlm(fψ,ψ[:,indices],',')
close(fψ)

fζ = open("./out/zeta.csv","w")
head = prod("zeta$i"*(i == indices[end] ? "\n" : ", ")  for i in indices)
write(fζ, head)
writedlm(fζ,ζ[:,indices],',')
close(fζ)

### Extract posterior summary measures, treating inital 50% of the chain as burnin;
# writing these measures to csv file
p = 0.05
A = view(ψ, size(ψ, 1)÷2:size(ψ, 1), :)
upper = mapslices(v-> quantile(v, 1 - p/2), A; dims= 1)
med = median(A, dims=1)
ave = mean(A,dims=1)
lower = mapslices(v-> quantile(v,  p/2), A; dims= 1)

summaryψ=[ave' lower' med' upper']
f = open("./out/summarypsi.csv","w")
write(f, "average, lower, median, upper \n")
writedlm(f,summaryψ,',')
close(f)

# also with p=0.5
p = 0.5
A = view(ψ, size(ψ, 1)÷2:size(ψ, 1), :)
upper = mapslices(v-> quantile(v, 1 - p/2), A; dims= 1)
med = median(A, dims=1)
ave = mean(A,dims=1)
lower = mapslices(v-> quantile(v,  p/2), A; dims= 1)

summaryψ=[ave' lower' med' upper']
f = open("./out/summarypsi50.csv","w")
write(f, "average, lower, median, upper \n")
writedlm(f,summaryψ,',')
close(f)
