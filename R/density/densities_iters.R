library(gridExtra)

data = as.numeric(read.csv("../../data.csv", sep=",", header=F))[1:100]
names(data)<-c("data")
density0.25 = read.csv("../../density0.25.csv", header=F)
dens.matrix = read.csv("../../dens_estimate_iterations.csv", header=F)

h<-hist(data, plot=F)
h$counts <- h$counts / sum(h$counts)

#x11()
pdf("densities_iters.pdf")
plot(h, freq=T, ylab="Densities", xlim=c(0,10), ylim=c(0,0.4),
  main="Local density estimates at single iterations")
for(i in 1:nrow(dens.matrix)){
  lines(density0.25[,1], dens.matrix[i,], col="gray", lwd = 1)
}
lines(density0.25[,1], density0.25[,2], col="black", lwd = 2)
legend("topright", c("Densities at iterations 1000,2000,...,15000",
  "Mean density"), lty=c(1,1), col=c("gray","black"))
dev.off()
q("no")
