generate.normal.samples <- function(mean1, mean2, sigma, numSamples) {
  r = runif(numSamples)
  component = as.integer(r > 0.5) + 1
  means = numeric(numSamples)
  means[component == 1] = mean1
  means[component == 2] = mean2
  return (rnorm(n=numSamples, mean=means, sd=sigma))
}

fit.stan.model <- function(H, num_samples, samples) {
  dat = list(
    H=H,
    max_num_samples = max(num_samples),
    num_samples=num_samples,
    samples=samples)

  model = stan_model(file = 'our_model.stan')
  fit <- sampling(model, data = dat, control = list(adapt_delta = 0.8), cores=4, iter=2000)
  return(fit)
}

visualize <- function(fit, samples, numSamples) {

  alphaPlot = plot(fit, plotfun = "trace", pars=c("alpha_0")) + ggtitle("Traceplot of alpha_0")

  predictions = extract(fit, "predictions")$predictions
  cols <- c("Group 1" = "red", "Group 2" = "green")
  hist = ggplot() +
    geom_histogram(mapping = aes(x = c(samples[1, 1:numSamples[1]], samples[2, 1:numSamples[2]]), y=..density..)) +
    labs(x="", y="") +
    geom_density(mapping = aes_(predictions[, 1], colour="Group 1")) +
    geom_density(mapping = aes_(predictions[, 2], colour="Group 2")) +
    scale_colour_manual(name="",values=cols) +
    geom_hline(yintercept=0, colour="white", size=0.5) +
    theme_bw() + ggtitle("Histogram of data and predictive densities")

  return(bayesplot_grid(alphaPlot, hist))
  # legend(x="topright", col=c("red", "forestgreen"), lty=c(2, 2), bty="n",
  #        legend=c("Posterior density group 1", "Posterior density group 2"))
}
