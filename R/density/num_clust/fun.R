library(ggplot2)

num_clust= as.numeric(read.csv("num_clust_best.csv", sep=",", header=F))[1:10]

mass= as.numeric(read.csv("values_mass.csv", sep=",", header=F))[1:10]
plot(mass,num_clust)

df=data.frame(mass,num_clust)

ggplot(df,aes(x=mass, y=num_clust))+
  geom_line(size=1)+
  geom_point(size=2)+
  ggtitle("Number of clusters as a function of the total mass") +
  theme_bw()+
  xlab("Total mass") + ylab("Number of clusters")
ggsave("num_clust.pdf")