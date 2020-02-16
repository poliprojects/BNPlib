library(ggplot2)
library(gridExtra)


data = as.numeric(read.csv("data.csv", sep=",", header=F))[1:100]
names(data)<-c("data")

density0.25 = read.csv("../../density0.25.csv", header=F) 

dens.matrix = read.csv("../../dens_estimate_iterations.csv", header=F)

h<-hist(data, plot=F)
h$counts <- h$counts / sum(h$counts)

#x11()
pdf("densities_iters.pdf")
plot(h, freq=TRUE, ylab="Densities", xlim=c(0,10), ylim=c(0,0.4),
  main="Predictive densities")
for(i in 1:nrow(dens.matrix)){
  lines(density0.25[,1], dens.matrix[i,], col="gray", lwd = 1)
}
lines(density0.25[,1], density0.25[,2], col="black", lwd = 2)
legend("topright", c("Densities at iterations 1000,2000,...,15000",
  "Mean density"), lty=c(1,1), col=c("gray","black"))
dev.off()


legend("bottomright",  c("Estimated posterior density with:", "M=0.35", "M=1.0", "M=5.0"), 
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


a=dnorm(density035[,1],mean=7.08222, sd=	sqrt(1.59500)) #0
b=dnorm(density035[,1],mean=	4.27202, sd=sqrt(1.40144)) #1
c=dnorm(density035[,1],mean=5.08181, sd=sqrt(2.04721)) #2
d=dnorm(density035[,1],mean=	2.19518, sd=sqrt(2.28840))#3


dat=data.frame(table(clust_final035[,3]))
dat[,2]=dat[,2]/100
#v=c(1:100)
#v[1:50]="red"
#v[51:100]="orange"

data=data.frame(data)
components=data.frame(density035[,1], a,b,c,d)
names(components)[names(components) == "density035...1."] <- "V1"

weight_comp=data.frame(table(clust_final035[,3])[1]*a, table(clust_final035[,3])[2]*b,table(clust_final035[,3])[3]*c,table(clust_final035[,3])[4]*d)/(100+0.35)
weight_comp=data.frame(density035[,1],weight_comp)
names(weight_comp)<-c("V1", "a", "b", "c", "d")

  
p1=ggplot(data, aes(x=data)) +geom_histogram(aes(y=..density..),binwidth=0.72, fill="white", color="black")+

geom_line(data = density1, aes(x = V1, y = V2,col = '5'), size=1)+
geom_line(data = weight_comp, aes(x = V1, y = a, col = '1')) +
geom_line(data = weight_comp, aes(x = V1, y = b, col = '2')) +
geom_line(data = weight_comp, aes(x = V1, y = c, col = '3')) +
geom_line(data = weight_comp, aes(x = V1, y = d, col = '4')) +

  ggtitle("Components of estimated posterior density")  +
  xlab("Data") +
  xlim(c(1,10))+
  ylab("Densities") +
  scale_color_manual(name="",values =c('1'="yellow",'2'= "orange", '3'= "green", '4'="blue", '5'="red"), 
                                     labels = c( "Density of cluster 0", "Density of cluster 1", "Density of cluster 2", "Density of cluster 3", "Estimated posterior density"))



p2=ggplot(data = dat, aes(x = reorder(Var1, -Freq), Freq)) + 
  geom_bar(stat = "identity", color='black',fill='gray') +
  ggtitle("Component weights") +
  xlab("Clusters") + ylab("Relative frequency")

#svg("components.svg")

grid.arrange(p1, p2, widths = c(5,2.5))
#dev.off()
