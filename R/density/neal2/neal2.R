library(ggplot2)

setwd("C:/Users/EzioGreggio/Desktop/density_prog/neal2")

data= data.frame(as.numeric(read.csv("data.csv", sep=",", header=F))[1:100])
densityneal2= read.csv("densityneal2.csv", sep=",", header=F) 
density0.25= read.csv("density0.25.csv", sep=",", header=F) 

names(data)<-c("data")




ggplot(data, aes(x=data)) +geom_histogram(aes(y=..density..),binwidth=0.72, fill="white", color="black")+
  geom_line(data = densityneal2, aes(x = V1, y = V2,col = '1'))+
  geom_line(data = density0.25, aes(x = V1, y = V2,col = '2'))+
  ggtitle("Estimated posterior density")+
  xlab("data") +
  xlim(c(1,10))+
  ylab("density") +scale_color_manual(name="",values =c('1'="blue",'2'="green"), 
                                      labels = c( "Neal2","Neal8"))


