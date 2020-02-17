library(ggplot2)
library(gridExtra)


data= data.frame(as.numeric(read.csv("data.csv", sep=",", header=F))[1:100])
names(data)<-c("data")


# iteration 2490 best with M=0.25 m=3
density_best_iter= read.csv("density_best_iter.csv", sep=",", header=F)
clust_best0.25= read.csv("clust_best0.25.csv", sep=",", header=F) 



#############################################
#m=3 M=0.25 with BEST clustering and density at best_clust iteration

a=dnorm(density_best_iter[,1],mean=4.00769, sd=	0.933806) #0
b=dnorm(density_best_iter[,1],mean=7.04811, sd=1.202110) #1



components=data.frame(density_best_iter[,1], a,b)
names(components)[names(components) == "density_best_iter...1."] <- "V1"

weight_comp=data.frame(table(clust_best0.25[,3])[1]*a, table(clust_best0.25[,3])[2]*b)/(100+0.25)
weight_comp=data.frame(density_best_iter[,1],weight_comp)
names(weight_comp)<-c("V1", "a", "b")

  
p1=ggplot(data, aes(x=data)) +geom_histogram(aes(y=..density..),binwidth=0.72, fill="white", color="black")+

geom_line(data = density_best_iter, aes(x = V1, y = V2,col = '3'), size=1)+
geom_line(data = weight_comp, aes(x = V1, y = a, col = '1')) +
geom_line(data = weight_comp, aes(x = V1, y = b, col = '2'))  +
  ggtitle("Components of the density at the best clustering")  +
  xlab("Data") +
  xlim(c(1,10))+
  ylab("Densities") +
  scale_color_manual(name="",values =c('1'="gray",'2'= "gray", '3'= "black"), 
                                     labels = c( "Density of cluster 0", "Density of cluster 1", "Full density"))




