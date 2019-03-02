using Test
using PointProcessInference

observations, parameters, extra = PointProcessInference.loadexample("coal")
res = PointProcessInference.inference(observations; parameters...)
PointProcessInference.showresults(res)
