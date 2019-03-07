
obs = rand(Exponential(2), 100)
res = PointProcessInference.inference(obs, summaryfile="info.txt")

PointProcessInference.ggplot2()
plotposterior(res)
