rm(list=ls())
library(rjags)
library(dclone)
library(foreach)
library(optparse)
library(parallel)
library(doMC)
# setwd("/home/mario/PhD/exchangeability/")
source("utils/data_utils.R")
source("utils/experiments_utils.R")

initialC = list()
initialC[[1]] = c(1, 2, 3, 3)
initialC[[2]] = c(4, 2, 3, 1)
initialC[[3]] = c(1, 1, 1, 1)
initialC[[4]] = c(1, 2, 3, 4)
initialC[[5]] = c(4, 4, 4, 4)

run.exp <- function(i) {
  if (i == 1) {
    ys = generate.example1(args$num_samples)
  } else if (i == 2) {
    ys = generate.example2(args$num_samples)
  } else if (i == 3) {
    ys = generate.example3(args$num_samples)
  } else if (i == 4) {
    ys = generate.example4(args$num_samples)
  } else if (i == 5){
    ys = generate.example5(args$num_samples)
  } else {
    stop("Experiment number must be 1, 2, 3 or 4")
  }
  
  print(sprintf("Running Experiment %i", i))
  outfile = file.path(outDir, sprintf("experiment%i_data.RData", i))
  saveRDS(ys, outfile)
  
  num_samples = rep(ncol(ys), nrow(ys))
  
  # one DP per group independently
  mcmc_samples = run.jags.parallel(
    "jags/semi_hdp.jags", nproc, ys, num_samples, H,
    inits=list(rho=initialC[[i]]))
  outfile = file.path(outDir, sprintf("experiment%i_semi_hdp.RData", i))
  saveRDS(mcmc_samples, outfile)
}


option_list = list(
  make_option(c("--out_dir"), type="character", default="/home/mario/PhD/exchangeability/results/test"),
  make_option(c("--num_samples"), type="integer", default=90)
)

opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser);


outDir = file.path(args$out_dir, Sys.Date())
nproc = 4
H = 20

if (dir.exists(outDir)) {
  outDir_ = outDir
  i = 2
  while (dir.exists(outDir_)) {
    outDir_ = paste(outDir, sprintf("v%i", i), sep = "_")
    i = i+1
  }
  outDir = outDir_
}

dir.create(outDir)

experiment_params = data.frame(list(
  "nproc" = nproc,
  "num_components" = H,
  "num_samples" = args$num_samples
))

write.csv(experiment_params, file = file.path(outDir, "parameters.csv"))


experiments = c(1, 2, 3, 4, 5)
experiments = c(3, 4, 5)

numCores = min(detectCores() - 1, length(experiments) * nproc) %/% nproc
print(sprintf("Registering %i cores at upper level, total number of used cores: %i",
              numCores, numCores * nproc))
cl <- makeCluster(numCores, type="FORK")

parLapply(cl, experiments, run.exp)

mcmc_samples = run.jags(
  "jags/univariate_eppf.txt", ys, num_samples, H, 1, 0, 100)
MCMCchains(mcmc_samples, params="maxRho")
