library(ggplot2)
library(gridExtra)

setwd("C:/Users/EzioGreggio/Desktop/density_prog/dens_withMm3")

data= data.frame(as.numeric(read.csv("data.csv", sep=",", header=F))[1:100])
names(data)<-c("data")

#varying total mass M

density0.25= read.csv("density0.25.csv", sep=",", header=F)
density5= read.csv("density5.csv", sep=",", header=F) 
density1= read.csv("density1.csv", sep=",", header=F) 
density0.5= read.csv("density0.5.csv", sep=",", header=F) 
density10= read.csv("density10.csv", sep=",", header=F) 



p <- ggplot(data, aes(x=data)) + 
  geom_density()

ggplot(data, aes(x=data)) +geom_histogram(aes(y=..density..),binwidth=0.72, fill="white", color="black")+
  geom_line(data = density0.25, aes(x = V1, y = V2,col = '1'))+
  geom_line(data = density0.5, aes(x = V1, y = V2,col = '2'))+
  geom_line(data = density1, aes(x = V1, y = V2,col = '3'))+
  geom_line(data = density5, aes(x = V1, y = V2,col = '4'))+
  geom_line(data = density10, aes(x = V1, y = V2,col = '5'))+
  
  ggtitle("Estimated posterior density")  +
  xlab("data") +
  xlim(c(1,10))+
  ylab("density") +scale_color_manual(name="",values =c('1'="red",'2'="green",'3'= "blue", '4'='orange', '5'='violet'), 
                                      labels = c( "M=0.25","M=0.5", "M=1.0", "M=5.0","M=10.0"))

