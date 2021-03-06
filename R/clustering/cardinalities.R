cards = read.csv("../../clust_cardinalities.csv", header=F)

# Pre-processing
cards = t( cards[1:(length(cards)-1)] )

#pdf("cardinalities_all.pdf")
#plot(cards)
#dev.off()

thin = seq(1, length(cards), 100)
pdf("cardinalities_thinned.pdf")
plot(thin, cards[thin],
	xlab="Iterations 100,200,...,15000",
	ylab="Number of clusters",
	main="Number of clusters at single iterations")
dev.off()
