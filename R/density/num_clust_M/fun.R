library(ggplot2)

num_clust= as.numeric(read.csv("num_clust_best.csv", sep=",", header=F))

mass= as.numeric(read.csv("mass_values.csv", sep=",", header=F))
plot(mass,num_clust)

df=data.frame(mass,num_clust)

ggplot(df,aes(x=mass, y=num_clust))+
  scale_x_log10()+
  geom_line(size=1)+
  geom_point(size=2)+
  ggtitle("Number of clusters as a function of the total mass") +
  theme_bw()+
  xlab("Total mass") + ylab("Number of clusters")
ggsave("num_clust.pdf")

