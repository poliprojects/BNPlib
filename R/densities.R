data  = read.csv("../data.csv", header=F)
dens  = read.csv("../density.csv", header=F)
cards = read.csv("../clust_cardinalities.csv", header=F)

# Pre-processing
data  = as.numeric(data[1:(length(data)-1)])
cards = cards[1:(length(cards)-1)]

#print(dim(cards))

svg(filename="density_pure.svg")
hist(data, freq=T)
lines(dens)
dev.off()

# TODO check if it works
#svg(filename="cardinalities_all.svg")
#plot(cards)
#dev.off()

#thin = seq(1,length(cards),500)
#print(thin)
#svg(filename="cardinalities_thinned.svg")
#plot(cards[thin])
#dev.off()
