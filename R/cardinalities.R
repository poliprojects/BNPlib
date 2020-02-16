data  = read.csv("../data.csv", header=F)
dens  = read.csv("../density0.25.csv", header=F)
cards = read.csv("../clust_cardinalities.csv", header=F)

# Pre-processing
data  = as.numeric(data[1:(length(data)-1)])
cards = t( cards[1:(length(cards)-1)] )

#pdf("cardinalities_all.pdf")
#plot(cards)
#dev.off()

thin = seq(1, length(cards), 100)
pdf("cardinalities_thinned.pdf")
plot(thin,cards[thin])
dev.off()
