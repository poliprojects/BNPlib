

data= as.numeric(read.csv("data.csv", sep=",", header=F))

density5= read.csv("density5.csv", sep=",", header=F) #24 clusters
density1= read.csv("density1.csv", sep=",", header=F) #9 clusters
density035= read.csv("density035.csv", sep=",", header=F) # 1 cluster


h<-hist(data, plot=F)
h$counts <- h$counts / sum(h$counts)

x11()
plot(h, freq=TRUE, ylab="Relative Frequency", ylim=c(0,0.4), xlim=c(0,10))

lines(density5[,1], density5[,2], col="green",  lwd = 1.5)
lines(density1[,1], density1[,2], col="red",  lwd = 1.5)
lines(density035[,1], density035[,2], col="blue",  lwd = 1.5)


legend("topright",  c("Estimated posterior density with:", "M=0.35", "M=1.0", "M=5.0"), 
       lty=c(0,1,1,1), 
       col=c("blue","red", "green"))




clust_best035= read.csv("clust_best035.csv", sep=",", header=F)
clust_final035= read.csv("clust_final035.csv", sep=",", header=F)

hist(clust_final035[,3])

clust_best1= read.csv("clust_best1.csv", sep=",", header=F) 
clust_final1= read.csv("clust_final1.csv", sep=",", header=F) 
hist(clust_best1[,3])
hist(clust_final1[,3])


library(ggplot2)

dat=data.frame(table(clust_best1[,3]))

ggplot(data = dat, aes(x = reorder(Var1, -Freq), Freq)) + 
  geom_bar(stat = "identity", color='skyblue',fill='steelblue') 


clust_best5= read.csv("clust_best5.csv", sep=",", header=F) 
clust_final5= read.csv("clust_final5.csv", sep=",", header=F) 



