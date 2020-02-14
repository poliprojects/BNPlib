data  = read.csv("../data.csv", header=F)
dens  = read.csv("../density.csv", header=F)
cards = read.csv("../clust_cardinalities.csv", header=F)

# Pre-processing
data  = as.numeric(data[1:(length(data)-1)])
cards = t( cards[1:(length(cards)-1)] )

svg(filename="density_pure.svg")
hist(data, freq=F, ylim=c(0.0, 0.3))
lines(dens)
#abline(v=4, col=1)
#abline(v=6, col=2)
dev.off()

svg(filename="cardinalities_all.svg")
plot(cards)
dev.off()

thin = seq(1, length(cards), 100)
svg(filename="cardinalities_thinned.svg")
plot(thin,cards[thin])
dev.off()
