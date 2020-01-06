rm(list=ls())
setwd("/home/mario/PhD/exchangeability/")
source("utils/univ_normal_utils.R")

library(sn)
library(rstan)
rstan_options(auto_write = TRUE)
rstan_options(mc.cores = parallel::detectCores() -2)
library(bayesplot)
library(ggplot2)

###########################
# Different Distributions #
###########################

num_samples = c(100, 100)
max_num_samples = max(num_samples)

samples1 = numeric(max_num_samples)
samples1[1:num_samples[1]] = generate.normal.samples(-10, -8, 0.5, num_samples[1])
# samples1[1:num_samples[1]] = rnorm(num_samples[1], -5, 0.5)

samples2 = numeric(max_num_samples)
samples2[1:num_samples[2]] = generate.normal.samples(12, 15, 0.5, num_samples[2])

samplesNoCommon = rbind(samples1, samples2)

ks.test(samples1, samples2)

fitNoCommon = fit.stan.model(5, num_samples, samplesNoCommon)  
visualize(fitNoCommon, samplesNoCommon, num_samples)

plot(fitNoCommon, plotfun = "hist", pars=c("isEx"))


#####################
# Same Distribution #
#####################

num_samples = c(500, 500)
max_num_samples = max(num_samples)

samples1 = numeric(max_num_samples)
samples1[1:num_samples[1]] = rnorm(num_samples[1], 0, 0.5)

samples2 = numeric(max_num_samples)
samples2[1:num_samples[2]] =rnorm(num_samples[2], 0, 0.5)

par(mfrow=c(2, 1))
hist(samples1)
hist(samples2)

samplesSame = rbind(samples1, samples2)

ks.test(samples1[1:num_samples[1]], samples2[1:num_samples[2]])

fitSame = fit.stan.model(5, num_samples, samplesSame)  
visualize(fitSame, samplesSame, num_samples)

plot(fitSame, plotfun = "hist", pars=c("isEx"))
