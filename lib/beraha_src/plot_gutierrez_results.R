############
# Analysis #
############

# rm(list=ls())
set.seed(123456)
library(bayesplot)

setwd("/home/mario/PhD/exchangeability/")
source("utils/data_utils.R")
source("utils/experiments_utils.R")


library(MCMCvis)
library(reticulate)
source_python("utils/python_utils.py")
source_python("utils/python_plots.py")
source("utils/plots.R")


plot.predictive.distrib <- function(mcmc_samples, observed_ys, experiment, expType, outfile) {
  if (expType == "shared") {
      y_ppd = get.predictive.samples(mcmc_samples, T)
    title = sprintf("Experiment %i - Shared atoms", experiment)
  } else if (expType == "base") {
    y_ppd = get.predictive.samples(mcmc_samples, F)
    title = sprintf("Experiment %i - Base", experiment)
  } else if (expType %in% c("mixing_var", "escobar_west", "semi_hdp")) {
    y_ppd = get.predictive.samples.mixing(mcmc_samples)
    title = sprintf("Experiment %i - Mixing Variance", experiment)
  } else
    stop("wrong experiment type")
  if (experiment == 1)
    plot_exp_1(y_ppd, observed_ys, title, outfile)
  else if (experiment == 2)
    plot_exp_2(y_ppd, observed_ys, title, outfile)
  else if (experiment == 3)
    plot_exp_3(y_ppd, observed_ys, title, outfile)
  else if (experiment == 4)
    plot_exp_4(y_ppd, observed_ys, title, outfile)
  else if (experiment == 5)
    plot_exp_5(y_ppd, observed_ys, title, outfile)

  return (y_ppd)
}


load.mcmc.samples <- function(experiment, expType, dir="results/gutierrez/") {
  if (!(expType %in% c("shared", "base", "mixing_var", "escobar_west", "semi_hdp")))
    stop("wrong experiment type")
  filename = paste(dir, sprintf("experiment%i_%s.RData", experiment, expType), sep="")
  if (!file.exists(filename))
    stop("file does not exist")
  return (readRDS(filename))
}

load.observed.data <- function(experiment, dir="results/gutierrez/") {
  filename = file.path(dir, sprintf("experiment%i_data.RData", experiment))
  return (readRDS(filename))
}

plot.all <- function(results_dir="results/gutierrez/", output_dir="report/images/") {
  correctClusters = list()
  correctClusters[[1]] = "[[1, 2, 3], [4]]"
  correctClusters[[2]] = "[[1, 4], [2], [3]]"
  correctClusters[[3]] = "[[1], [2], [3], [4]]"
  correctClusters[[4]] = "[[1, 2, 3, 4]]"
  correctClusters[[5]] = "[[1], [2], [3], [4]]"
  expDate = strsplit(results_dir, "/")[[1]][3]
  for (experiment in c(1:5)) {
    correctClus = correctClusters[[experiment]]
    for (expType in c("shared", "base", "mixing_var", "escobar_west")) {
      tryCatch({
        mcmc_samples = load.mcmc.samples(experiment, expType, results_dir)
        observed = load.observed.data(experiment, results_dir)
        outfile = paste(output_dir,
                        sprintf("gutierrez_exp%i_%s_clust_%s.png", experiment, expType, expDate),
                        sep="")
        barplot.clusters(MCMCchains(mcmc_samples, params="rho"),
                         correctClusters[[experiment]],
                         outfile)
        
        outfile = paste(output_dir,
                        sprintf("gutierrez_exp%i_%s_gammachain_%s.png", experiment, expType, expDate),
                        sep="")
        gammaChain = MCMCchains(mcmc_samples, params="gam")
        png(outfile)
        mcmc_trace(gammaChain)
        dev.off()

        outfile = paste(output_dir,
                        sprintf("gutierrez_exp%i_%s_predictive_%s.png", experiment, expType, expDate),
                        sep="")
        plot.predictive.distrib(mcmc_samples, observed, experiment, expType, outfile)
      }, error = function(err_cond) {
        cat(sprintf(
          "Could not process experiment %i, type: %s, probably the file is not present\n",
          experiment, expType))
      })
    }
  }
}

plot.all("results/gutierrez/2019-04-03_v2/")



correctClusters = list()
correctClusters[[1]] = "[[1, 2, 3], [4]]"
correctClusters[[2]] = "[[1, 4], [2], [3]]"
correctClusters[[3]] = "[[1], [2], [3], [4]]"
correctClusters[[4]] = "[[1, 2, 3, 4]]"
correctClusters[[5]] = "[[1], [2], [3], [4]]"

results_dir = "results/gutierrez/2019-05-10/"
expType = "semi_hdp"


experiment = 3
mcmc_samples = load.mcmc.samples(experiment, expType, results_dir)
observed = load.observed.data(experiment, results_dir)
barplot.clusters(MCMCchains(mcmc_samples, params="rho"),
                 correctClusters[[experiment]])

gammaChain = MCMCchains(mcmc_samples, params="gam")
mcmc_trace(gammaChain)

a = plot.predictive.distrib(mcmc_samples, observed, experiment, expType, "")
