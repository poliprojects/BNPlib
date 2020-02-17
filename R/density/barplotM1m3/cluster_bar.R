library(ggplot2)
library(gridExtra)


data= data.frame(as.numeric(read.csv("data.csv", sep=",", header=F))[1:100])
names(data)<-c("data")


#############################################
#BARPLOT
#M=1 with best clustering

clust_best1= read.csv("clust_best1.csv", sep=",", header=F) 

dat=data.frame(table(clust_best1[,3]))
dat[,2]=dat[,2]/100

p2=ggplot(data = dat, aes(x = reorder(Var1, -Freq), Freq)) + 
  geom_bar(stat = "identity", color='black',fill='gray') +
  ggtitle("Component weights") +
  xlab("Clusters") + ylab("Relative frequency")


