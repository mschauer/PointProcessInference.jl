module PointProcessInference

using Distributions
using Optim
using Random, Statistics
using DelimitedFiles
using SpecialFunctions
using DataDeps
using DataFrames # FIX: add as dependency

lgamma(x) = (logabsgamma(x))[1]
include("funs-ig.jl")
include("gen-extract-data.jl")
include("revjump.jl")

export mloglikelihood, computebinning, revjump,  evalstepfunction

function __init__()
    registerdatadeps()
end

include("ppp.jl")

plotscript() = joinpath(dirname(pathof(PointProcessInference)), "..", "contrib", "process-output-simple.jl")

end # module
