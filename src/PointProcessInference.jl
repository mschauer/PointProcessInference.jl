module PointProcessInference
using DataFrames
using Distributions
using Optim
using Random, Distributions, Statistics
using DelimitedFiles
using SpecialFunctions
using RCall
using DataFrames
using TimerOutputs
using DataDeps

include("funs-ig.jl")
include("process-output-simple.jl")
include("gen-extract-data.jl")

function __init__()
    registerdatadeps()
end

include("ppp.jl")

end # module
