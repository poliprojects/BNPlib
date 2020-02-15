

data= as.numeric(read.csv("data.csv", sep=",", header=F))

density5= read.csv("density5.csv", sep=",", header=F) #24 clusters
density1= read.csv("density1.csv", sep=",", header=F) #9 clusters
density045= read.csv("density045.csv", sep=",", header=F) # 3 clusters
density035= read.csv("density035.csv", sep=",", header=F) # 1 cluster


h<-hist(data, plot=F)
h$counts <- h$counts / sum(h$counts)

x11()
plot(h, freq=TRUE, ylab="Relative Frequency", ylim=c(0,0.4), xlim=c(0,10))

lines(density5[,1], density5[,2], col="green",  lwd = 1.5)
lines(density1[,1], density1[,2], col="red",  lwd = 1.5)
#lines(density045[,1], density045[,2], col="green",  lwd = 1.5)
lines(density035[,1], density035[,2], col="blue",  lwd = 1.5)


legend("topright",  c("Estimated posterior density with:", "M=0.35", "M=1.0", "M=5.0"), 
       lty=c(0,1,1,1), 
       col=c("blue","red", "green"))





