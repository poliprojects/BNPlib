library(ggplot2)


data= data.frame(as.numeric(read.csv("data.csv", sep=",", header=F))[1:100])
densityneal2= read.csv("densityneal2.csv", sep=",", header=F) 
density0.25= read.csv("density0.25.csv", sep=",", header=F) 

names(data)<-c("data")

x11()
plot(h, freq=TRUE, ylab="Relative Frequency", ylim=c(0,0.4), xlim=c(0,10))
lines(densityneal2[,1], densityneal2[,2], col="green",  lwd = 1.5)
lines(density0.25[,1], density0.25[,2], col="red",  lwd = 1.5)




ggplot(data, aes(x=data)) +geom_histogram(aes(y=..density..),binwidth=0.72, fill="white", color="black")+
  geom_line(data = densityneal2, aes(x = V1, y = V2,col = '1'))+
  geom_line(data = density0.25, aes(x = V1, y = V2,col = '2'))+
  ggtitle("Estimated posterior density")+
  xlab("data") +
  ylab("density") +scale_color_manual(name="",values =c('1'="red",'2'="green"), 
                                      labels = c( "Neal2","Neal8"))


