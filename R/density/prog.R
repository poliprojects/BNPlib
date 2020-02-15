library(ggplot2)
library(gridExtra)


data= data.frame(as.numeric(read.csv("data.csv", sep=",", header=F))[1:100])
names(data)<-c("data")

density5= read.csv("density5.csv", sep=",", header=F) 
density1= read.csv("density1.csv", sep=",", header=F) 
density05= read.csv("density05.csv", sep=",", header=F) 
density035= read.csv("density035.csv", sep=",", header=F) 


h<-hist(data, plot=F)
h$counts <- h$counts / sum(h$counts)

x11()
plot(h, freq=TRUE, ylab="Relative Frequency", ylim=c(0,0.4), xlim=c(0,10))
lines(density035[,1], density035[,2], col="green",  lwd = 1.5)
#lines(density05[,1], density05[,2], col="blue",  lwd = 1.5)
lines(density1[,1], density1[,2], col="orange",  lwd = 1.5)
lines(density5[,1], density5[,2], col="red",  lwd = 1.5)


legend("topright",  c("Estimated posterior density with:", "M=0.35", "M=1.0", "M=5.0"), 
       lty=c(0,1,1,1), 
       col=c("green","orange", "red"))



ggplot(data, aes(x=data)) +geom_histogram(aes(y=..density..),binwidth=0.72, fill="white", color="black")+
  
geom_line(data = density035, aes(x = V1, y = V2,col = '1'))+
geom_line(data = density1, aes(x = V1, y = V2,col = '2'))+
geom_line(data = density5, aes(x = V1, y = V2,col = '3'))+
ggtitle("Estimated posterior density")  +
xlab("data") +
ylab("density") +scale_color_manual(name="",values =c('1'="red",'2'="green",'3'= "blue"), 
                                                 labels = c( "M=0.35", "M=1.0", "M=5.0"))


              

clust_best035= read.csv("clust_best035.csv", sep=",", header=F)
clust_final035= read.csv("clust_final035.csv", sep=",", header=F)


clust_best1= read.csv("clust_best1.csv", sep=",", header=F) 
clust_final1= read.csv("clust_final1.csv", sep=",", header=F) 



clust_best05= read.csv("clust_best05.csv", sep=",", header=F)
clust_final05= read.csv("clust_final05.csv", sep=",", header=F)


clust_best5= read.csv("clust_best5.csv", sep=",", header=F) 
clust_final5= read.csv("clust_final5.csv", sep=",", header=F) 

#############################################
#M=1 with final


a=dnorm(density1[,1],mean=	4.81909, sd=	1.672570)#0
b=dnorm(density1[,1],mean=	5.26594, sd=	2.783170) #1
c=dnorm(density1[,1],mean=5.17475, sd	=1.166610) #2
d=dnorm(density1[,1],mean=	4.67092, sd=	0.349494)#3


dat=data.frame(table(clust_final1[,3]))
dat[,2]=dat[,2]/100
#v=c(1:100)
#v[1:50]="red"
#v[51:100]="orange"

data=data.frame(data)
components=data.frame(density1[,1], a,b,c,d)
names(components)[names(components) == "density1...1."] <- "V1"



p1=ggplot(data, aes(x=data)) +geom_histogram(aes(y=..density..),binwidth=0.75, fill="white", color="black")+

geom_line(data = density1, aes(x = V1, y = V2,col = '1'))+
geom_line(data = components, aes(x = V1, y = a, col = '2')) +
geom_line(data = components, aes(x = V1, y = b, col = '3')) +
geom_line(data = components, aes(x = V1, y = c, col = '4')) +
geom_line(data = components, aes(x = V1, y = d, col = '5')) +

  ggtitle("Estimated posterior density")  +
  xlab("data") +
  ylab("Relative frequency") +scale_color_manual(name="",values =c('1'="red",'2'="yellow",'3'= "orange", '4'= "green", '5'="blue"), 
                                     labels = c("Estimated posterior density", "Density of cluster 0", "Density of cluster 1", "Density of cluster 2", "Density of cluster 3"))


p2=ggplot(data = dat, aes(x = reorder(Var1, -Freq), Freq)) + 
  geom_bar(stat = "identity", color='black',fill='gray') +
  ggtitle("Barplot of Clusters") +
  xlab("Clusters") + ylab("Relative frequency")


grid.arrange(p1, p2, widths = c(3,2))