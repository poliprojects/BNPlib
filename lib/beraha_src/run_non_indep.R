rm(list=ls())
library(rjags)
library(dclone)
library(foreach)
library(optparse)
library(parallel)
library(doMC)
library(pbapply)
setwd("/home/mario/PhD/exchangeability/")
source("utils/data_utils.R")
source("utils/experiments_utils.R")

library(devtools)
# compileAttributes(pkgdir = ".", verbose = getOption("verbose"))
# pkgbuild::clean_dll()
devtools::load_all()

run.exp <- function(dat_) {
  i = dat_$i
  ys = dat_$ys
  print(sprintf("Running Experiment %i", i))
  outfile = file.path(outDir, sprintf("experiment%i_data.RData", i))
  saveRDS(ys, outfile)
  
  num_samples = rep(ncol(ys), nrow(ys))
  
  mcmc_samples = run.jags.parallel(
    "jags/semi_hdp.jags", nproc, ys, num_samples, H, 
    inits=list(rho=c(1, 2)))
  outfile = file.path(outDir, sprintf("experiment%i_semi_hdp.RData", i))
  saveRDS(mcmc_samples, outfile)
}

run.exp.gibbs <- function(dat_) {
  i = dat_$i
  ys = dat_$ys
  data = list()
  data[[1]] = ys[1, ]
  data[[2]] = ys[2, ]
  
  print(sprintf("Running Experiment %i", i))
  outfile = file.path(outDir, sprintf("experiment%i_data.RData", i))
  saveRDS(ys, outfile)
  
  num_samples = rep(ncol(ys), nrow(ys))
  
  outfile = file.path(outDir, sprintf("semi_hdp_experiment%i_chains.RData", i))
  chains = runGibbs("marginal_semi_hdp", data, 10000, 10000, 10000, 10, 10, 500, T)
  saveRDS(chains, outfile)

  # outfile = file.path(outDir, sprintf("semi_hdp_experiment%i_chains.RData", i))
  # chains = runGibbs("semi_hdp", data, 10000, 10000, 10000, 10, 10, 500, T)
  # saveRDS(chains, outfile)
  # 
  # outfile = file.path(outDir, sprintf("normal_conjugate_experiment%i_chains.RData", i))
  # chains = runGibbs("normal_conjugate", data, 10000, 10000, 10000, 10, 10, 500, T)
  # saveRDS(chains, outfile)
}


option_list = list(
  make_option(c("--out_dir"), type="character", default="/home/mario/PhD/exchangeability/results/non_indep"),
  make_option(c("--numcores"), type="integer", default=7)
)

opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser);


outDir = file.path(args$out_dir, Sys.Date())

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


experiments = 1:50
dat_ = list()
for (exp in experiments) {
  elem = list()
  elem$i = exp
  elem$ys = generate.non.indep3(100)
  dat_[[exp]] = elem
}

numCores = min(args$numcores, length(experiments) * nproc) %/% nproc
print(sprintf("Registering %i cores at upper level, total number of used cores: %i",
              numCores, numCores * nproc))

# cl <- makeCluster(numCores, type="FORK")
cl <- makeCluster(args$numcores, type="FORK")
pblapply(dat_, run.exp.gibbs, cl=cl)
