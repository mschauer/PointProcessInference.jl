module PointProcessInference

using DataFrames
using Distributions
using Optim
using Random
using DelimitedFiles
using SpecialFunctions
using RCall
using DataFrames
using TimerOutputs

include("funs-ig.jl")
include("process-output-simple.jl")
include("ppp.jl")

end # module
