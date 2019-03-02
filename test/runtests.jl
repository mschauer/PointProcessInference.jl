using Test
using PointProcessInference

observations, parameters, λinfo = PointProcessInference.loadexample("generated")
res = PointProcessInference.inference(observations; parameters...)
#include(joinpath(dirname(pathof(PointProcessInference)), "..", "contrib", "process-output-simple.jl")
