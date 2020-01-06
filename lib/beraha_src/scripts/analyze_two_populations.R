rm(list=ls())

setwd("/home/mario/PhD/exchangeability/scripts/")

analyzeExp <- function(expnum, outdir) {
  c = read.csv(sprintf(outdir + "exp%i_c_blocked.csv", expnum), header=F)
  preds = read.csv(sprintf(outdir + "exp%i_preds_blocked.csv", expnum), header=F)
}

outdir = "/home/mario/PhD/exchangeability/results/two_pop/"
expnum = 3
# cBlocked = read.csv(paste(outdir, sprintf("blocked_exp%i_c.csv", expnum), sep=""), header=F)
# predsBlocked = read.csv(paste(outdir, sprintf("blocked_exp%i_preds.csv", expnum), sep=""), header=F)
cMarginal = read.csv(paste(outdir, sprintf("marginal_exp%i_c.csv", expnum), sep=""), header=F)
predsMarginal = read.csv(paste(outdir, sprintf("marginal_exp%i_preds.csv", expnum), sep=""), header=F)

# par(mfrow=c(1, 2))
# plot(density(predsBlocked[, 1]))
# plot(density(predsBlocked[, 2]))

par(mfrow=c(1, 2))
plot(density(predsMarginal[, 1]))
plot(density(predsMarginal[, 2]))

# mean(cBlocked[, 1] == cBlocked[, 2])
mean(cMarginal[, 1] == cMarginal[, 2])



#########
# META  #
#########
library(ggplot2)

probas = numeric(50)
for (i in 1:50) {
  cChain = read.csv(sprintf("/home/mario/PhD/exchangeability/results/blocked/non_indep/exp%i_c_blocked.csv", i-1), header=F)
  probas[i] = mean(cChain[, 1] == cChain[, 2])
}
probas
hist(probas)
