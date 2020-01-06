setwd("/home/mario/PhD/exchangeability/")

library(sn)
library(rjags)
library(coda)
library(bayesplot)
library(ggplot2)
library(MCMCvis)
library(dclone)

source("utils/data_utils.R")
source("utils/experiments_utils.R")

library(MCMCvis)
library(reticulate)
source_python("utils/python_utils.py")
source_python("utils/python_plots.py")
source("utils/plots.R")

H = 20

###########################
# Different Distributions #
###########################

num_samples = c(100, 100)
max_num_samples = max(num_samples)

samples1 = numeric(max_num_samples)
samples1[1:num_samples[1]] = normal.mixture.samples(c(-10, -8), c(0.5, 0.5), 0.5, num_samples[1])
# samples1[1:num_samples[1]] = rnorm(num_samples[1], -5, 0.5)

samples2 = numeric(max_num_samples)
samples2[1:num_samples[2]] = normal.mixture.samples(c(12, 15), c(0.5, 0.5), 0.5, num_samples[2])
ys = rbind(samples1, samples2)

ks.test(samples1, samples2)

samplesDifferentAtoms = run.jags("jags/base_univariate.txt", ys, num_samples, H)

# yppdChain = MCMCchains(samplesDifferentAtoms, params="y_ppd")
rhoChain = MCMCchains(samplesDifferentAtoms, params="rho")
sigmaHChain = MCMCchains(samplesDifferentAtoms, params="sigmaH")
mcmc_trace(sigmaHChain)

is_equal = as.integer(rhoChain[, 1] == rhoChain[, 2])
is_equal_tab = table(is_equal)
plot(prop.table(is_equal_tab), lwd=10)
title(main = "Posterior distribution of (c_1 = c_2)")
mtext("Different distributions")


samplesSharedAtoms = run.jags("jags/univariate_shared_atoms.txt", ys, num_samples, H)

yppdChain = MCMCchains(samplesSharedAtoms, params="y_ppd")
rhoChain = MCMCchains(samplesSharedAtoms, params="rho")
sigmaHChain = MCMCchains(samplesDifferentAtoms, params="sigmaH")
is_equal = as.integer(rhoChain[, 1] == rhoChain[, 2])
is_equal_tab = table(is_equal)
plot(prop.table(is_equal_tab), lwd=10)
title(main = "Posterior distribution of (c_1 = c_2) - Shared atoms")
mtext("Different distributions")

#####################
# Same Distribution #
#####################

num_samples = c(100, 100)
max_num_samples = max(num_samples)

samples1 = numeric(max_num_samples)
samples1[1:num_samples[1]] = rnorm(num_samples[1], 0, 0.5)

samples2 = numeric(max_num_samples)
samples2[1:num_samples[2]] =rnorm(num_samples[2], 0, 0.5)

ys = rbind(samples1, samples2)

ks.test(samples1[1:num_samples[1]], samples2[1:num_samples[2]])

samplesDifferentAtoms = run.jags("jags/base_univariate.txt", ys, num_samples, H)

yppdChain = MCMCchains(samplesDifferentAtoms, params="y_ppd")
rhoChain = MCMCchains(samplesDifferentAtoms, params="rho")
sigmaHChain = MCMCchains(samplesDifferentAtoms, params="sigmaH")
is_equal = as.integer(rhoChain[, 1] == rhoChain[, 2])
is_equal_tab = table(is_equal)
plot(prop.table(is_equal_tab), lwd=10)
title(main = "Posterior distribution of (c_1 = c_2)")
mtext("Equal distributions")

mcmc_trace(sigmaHChain)

samplesSharedAtoms = run.jags("jags/univariate_shared_atoms.txt", ys, num_samples, H)

yppdChain = MCMCchains(samplesSharedAtoms, params="y_ppd")
rhoChain = MCMCchains(samplesSharedAtoms, params="rho")
sigmaHChain = MCMCchains(samplesSharedAtoms, params="sigmaH")
mu0Chain = MCMCchains(samplesSharedAtoms, params="mu0")
is_equal = as.integer(rhoChain[, 1] == rhoChain[, 2])
is_equal_tab = table(is_equal)
plot(prop.table(is_equal_tab), lwd=10)
title(main = "Posterior distribution of (c_1 = c_2) - Shared atoms")
mtext("Equal distributions")

mcmc_trace(mu0Chain)

########################
# One common component #
########################

num_samples = c(100, 100)
max_num_samples = max(num_samples)

samples1 = numeric(max_num_samples)
samples1[1:num_samples[1]] = normal.mixture.samples(c(-5, 0), c(0.5, 0.5), 0.5, num_samples[1])


samples2 = numeric(max_num_samples)
samples2[1:num_samples[2]] = normal.mixture.samples(c(0, 10), c(0.5, 0.5), 0.5, num_samples[2])
ys = rbind(samples1, samples2)

par(mfrow=c(2,1))
hist(samples1, xlim=c(-6, 15), prob=TRUE, nclass = 20)
hist(samples2, xlim=c(-6, 15), prob=TRUE, nclass=20)
dev.off()

ks.test(samples1, samples2)

mcmc_samples = run.jags.parallel(
  "jags/univariate_like_escobar_west.jags", 4, ys, num_samples, H, inits=list("rho"=c(2, 1)))

barplot.clusters(MCMCchains(mcmc_samples, params="rho"), "[[1], [2]]")
