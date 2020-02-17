library(ggplot2)
library(gridExtra)

setwd("C:/Users/EzioGreggio/Desktop/density_prog/dens_with_mM0.25")

data= data.frame(as.numeric(read.csv("data.csv", sep=",", header=F))[1:100])
names(data)<-c("data")

#varying m aux variables 

density_m1= read.csv("density_m1.csv", sep=",", header=F)
density_m3= read.csv("density_m3.csv", sep=",", header=F) 
density_m5= read.csv("density_m5.csv", sep=",", header=F) 
density_m10= read.csv("density_m10.csv", sep=",", header=F) 
density_m50= read.csv("density_m50.csv", sep=",", header=F) 


p <- ggplot(data, aes(x=data)) + 
  geom_density()

ggplot(data, aes(x=data)) +geom_histogram(aes(y=..density..),binwidth=0.72, fill="white", color="black")+
  geom_line(data = density_m1, aes(x = V1, y = V2,col = '1'))+
  geom_line(data = density_m3, aes(x = V1, y = V2,col = '2'))+
  geom_line(data = density_m5, aes(x = V1, y = V2,col = '3'))+
  geom_line(data = density_m10, aes(x = V1, y = V2,col = '4'))+
  geom_line(data = density_m50, aes(x = V1, y = V2,col = '5'))+
  xlim(c(1,10))+
  ggtitle("Estimated posterior density")  +
  xlab("data") +
  ylab("density") +scale_color_manual(name="",values =c('1'="red",'2'="green",'3'= "blue", '4'='orange', '5'='violet'), 
                                      labels = c( "m=1","m=3", "m=5", "m=10","m=50"))


