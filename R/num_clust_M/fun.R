library(ggplot2)

num_clust_best= as.numeric(read.csv("num_clust_best.csv", sep=",", header=F))[1:10]
#num_clust_fin= as.numeric(read.csv("num_clusters_final.csv", sep=",", header=F))[1:10]

num_clust=num_clust_best
mass= as.numeric(read.csv("values_mass.csv", sep=",", header=F))[1:10]
plot(mass,num_clust)

df=data.frame(mass,num_clust)
ggplot(df,aes(x=mass, y=num_clust))+
  geom_line(linetype = "dashed")+
  geom_point()+
  ggtitle("Number of clusters in function of the total mass M") +
  xlab("M") + ylab("n. clusters")
