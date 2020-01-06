library(reticulate)
require(MCMCvis)
source("utils/experiments_utils.R")
source_python("utils/python_plots.py")


plot.predictive.distrib <- function(experiment, shared, results_dir, outfile) {
  y_ppd = get.predictive.samples(experiment, shared, results_dir)
  if (shared)
    title = sprintf("Experiment %i - Shared atoms", experiment)
  else
    title = sprintf("Experiment %i - Base", experiment)
  if (experiment == 1)
    plot_exp_1(y_ppd, title, outfile)
  else if (experiment == 2)
    plot_exp_2(y_ppd, title, outfile)
  else if (experiment == 3)
    plot_exp_3(y_ppd, title, outfile)
  else if (experiment == 4)
    plot_exp_4(y_ppd, title, outfile)
}
