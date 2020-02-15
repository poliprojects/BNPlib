library(ggplot2)


data= as.numeric(read.csv("data.csv", sep=",", header=F))[1:100]
density= read.csv("density.csv", sep=",", header=F) 

h<-hist(data, plot=F)
h$counts <- h$counts / sum(h$counts)

best_clust= read.csv("best_clust.csv", sep=",", header=F)
final_clust= read.csv("final_clust.csv", sep=",", header=F)

hist(clust_final035[,3])

x11()
plot(h, freq=TRUE, ylab="Relative Frequency", ylim=c(0,0.4), xlim=c(0,10))
lines(density[,1], density[,2], col="green",  lwd = 1.5)

