require(dclone)
library(mcclust)
library(parallel)
require(rjags)

run.jags.pararallel.helper <- function(i, data_, filename, inits=NULL) {
  var = c(
    "gam", "rho", "sigma", "means", "weights", "mu0", "sigmaH", "lambda",
    "probas", "wComp", "w", "GtildeAtoms", "weights00", "means00")
  if (is.null(inits)) {
    mod = jags.model(file=filename, data=data_, n.chains=1, n.adapt=5000)
  } else {
    mod = jags.model(file=filename, data=data_, n.chains=1, n.adapt=5000,
                     inits=inits)
  }
  samples = coda.samples(mod, var, 10000, thin = 5)
  return(samples[[1]])
}

run.jags.parallel <- function(filename, nproc, obs, num_samples, H, inits=NULL) {
  num_groups = nrow(obs)
  data_ = list(y=obs, H=H, num_samples=num_samples, numGroups=num_groups)
  set.seed(123456789, kind ="L'Ecuyer-CMRG");
  samples = mclapply(
      1:nproc, run.jags.pararallel.helper, mc.cores=nproc,
      data_=data_, filename=filename, inits=inits)
  return(do.call(rbind, samples))
}

run.jags <- function(filename, obs, num_samples, H, inits=NULL,
                     nadapt=5000, burnin=5000, nsamples=5000) {
   var = c(
     "gam", "rho", "sigma", "means", "weights", "mu0", "sigmaH", "lambda",
     "probas", "w00", "w", "GtildeAtoms", "weights00", "means00")
  num_groups = nrow(obs)
  data_ = list(y=obs, H=H, num_samples=num_samples, numGroups=num_groups)
  if (is.null(inits)) {
    mod = jags.model(file=filename, data=data_, n.chains=1, n.adapt=nadapt)
  } else {
    mod = jags.model(file=filename, data=data_, n.chains=1, n.adapt=nadapt,
                     inits=inits)
  }
  if (burnin > 0)
    update(mod, n.iter=burnin)
  samples = coda.samples(
    mod, var, nsamples, thin = 5)
}


load.mcmc.samples <- function(experiment, expType, dir="results/gutierrez/") {
  filename = paste(dir, sprintf("experiment%i_%s.RData", experiment, expType), sep="")
  return (readRDS(filename))
}

get.predictive.samples <- function(mcmc_samples, shared, numGroups=4, H=20) {
    rhochain = MCMCchains(mcmc_samples, params="rho")
    weightsChain = MCMCchains(mcmc_samples, params="weights")
    meansChain = MCMCchains(mcmc_samples, params="means")
    sigmaChain = MCMCchains(mcmc_samples, params="sigma")
    y_ppd = matrix(nrow=nrow(rhochain), ncol=numGroups)
    for (j in 1:numGroups) {
        for (i in 1:nrow(rhochain)) {
            rho = as.vector(rhochain[i, ])
            means = matrix(meansChain[i,], nrow=4)
            weights = matrix(weightsChain[i,], nrow=4)
            sigma = sigmaChain[i]
            s = numeric(numGroups)
            for (g in 1:numGroups) {
                s[g] = sample(c(1:H), 1, prob=weights[g, ])
            }
            if (shared)
              y_ppd[i, j] = rnorm(1, mean=means[s[rho[j]]], sigma)
            else
              y_ppd[i, j] = rnorm(1, mean=means[rho[j], s[rho[j]]], sigma)
        }
    }
    return(y_ppd)
}

binder.cluster <- function(experiment, shared, results_dir) {
    mcmc_samples = load.mcmc.samples(experiment, shared, results_dir)
    chain = MCMCchains(mcmc_samples, params="rho")
    clus = minbinder(comp.psm(chain), cls.draw=chain, method="all")$cl[1, ]
    return(indicatorToCluster(clus))
}

get.predictive.samples.mixing <- function(mcmc_samples, numGroups=4, H=20) {
    rhochain = MCMCchains(mcmc_samples, params="rho")
    weightsChain = MCMCchains(mcmc_samples, params="weights")
    meansChain = MCMCchains(mcmc_samples, params="means")
    sigmaChain = MCMCchains(mcmc_samples, params="sigma")
    y_ppd = matrix(nrow=nrow(rhochain), ncol=numGroups)
    for (j in 1:numGroups) {
        for (i in 1:nrow(rhochain)) {
            rho = as.vector(rhochain[i, ])
            means = matrix(meansChain[i,], nrow=numGroups)
            weights = matrix(weightsChain[i,], nrow=numGroups)
            stdevs = matrix(sigmaChain[i, ], nrow=numGroups)
            s = numeric(numGroups)
            for (g in 1:numGroups) {
                s[g] = sample(c(1:H), 1, prob=weights[g, ])
            }
            y_ppd[i, j] = rnorm(1, mean=means[rho[j], s[rho[j]]], stdevs[rho[j], s[rho[j]]])
        }
    }
    return(y_ppd)
}
