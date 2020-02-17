#data      = read.csv("../data.csv", header=F)
diss.mean = read.csv("../../dissim_matr_mean.csv", header=F)
diss.best = read.csv("../../dissim_matr_best.csv", sep=" ", header=F)

# Pre-processing
#data = as.numeric(data[1:(length(data)-1)])
diss.mean = diss.mean[,1:(length(diss.mean)-1)]
diss.mean = as.matrix(diss.mean)
diss.best = as.matrix(diss.best)

pdf("diss_mean.pdf")
image(diss.mean)
dev.off()

pdf("diss_best.pdf")
image(diss.best)
dev.off()

#diss.mean.round = round(diss.mean)
#pdf("diss_mean_round.pdf")
#image(diss.mean.round)
#dev.off()

#pdf("diss_comparison.pdf")
#image(abs(diss.mean.round-diss.best))
#dev.off()
