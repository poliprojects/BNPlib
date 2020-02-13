data = read.csv("density.csv", header=F)
plot(data)

#clust = read.csv("best_clust.csv")

x = seq(1, 6.0, 0.1)
data = dnorm(x, 4, 2/1.5)
#data = dinvgamma(x, 2.0, 2.0)
#plot(data)
