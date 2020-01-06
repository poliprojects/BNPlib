require(sn)

set.seed(20190402)

normal.mixture.samples <- function(means, weights, sigma, numSamples) {
  means_ = sample(x=means, numSamples, replace=T, prob=weights)
  return (rnorm(n=numSamples, mean=means_, sd=sigma))
}

skew.normal.samples <- function(mu, tau, alpha, numSamples) {
  return (sn::rmsn(n=numSamples, xi=mu, tau, alpha))
}

generate.example1 <- function(numSamples = 90) {
  out = matrix(nrow=4, ncol=numSamples)
  out[1, ] = skew.normal.samples(0, 1, 0, numSamples)
  out[2, ] = skew.normal.samples(0, 1, 0, numSamples)
  out[3, ] = skew.normal.samples(0, 1, 0, numSamples)
  out[4, ] = skew.normal.samples(0, 1, 1, numSamples)
  return(out)
}

generate.example2 <- function(numSamples = 90) {
  out = matrix(nrow=4, ncol=numSamples)
  out[1, ] = rnorm(numSamples, 0, 1)
  out[2, ] = rnorm(numSamples, 0, 2.25)
  out[3, ] = rnorm(numSamples, 0, 0.25)
  out[4, ] = rnorm(numSamples, 0, 1)
  return(out)
}

generate.example3 <- function(numSamples = 90) {
  out = matrix(nrow=4, ncol=numSamples)
  out[1, ] = rnorm(numSamples, 0, 0.49)
  out[2, ] = rnorm(numSamples, 0.6, 1)
  out[3, ] = normal.mixture.samples(c(-1.2, 1.2), c(0.5, 0.5), 0.25, numSamples)
  out[4, ] = normal.mixture.samples(c(-1.2, 1.2), c(0.3, 0.7), 0.25, numSamples)
  return(out)
}

# generate.example4 <- function(numSamples = 90) {
#   out = matrix(nrow=4, ncol=numSamples)
#   for (i in 1:4)
#     out[i, ] = normal.mixture.samples(c(-1.2, 1.2), c(0.5, 0.5), 0.25, numSamples)
#   return(out)
# }

generate.example4 <- function(numSamples = 90) {
  out = matrix(nrow=4, ncol=numSamples)
  for (i in 1:4)
    out[i, ] = rt(df=15, n=numSamples)
  return(out)
}


generate.example5 <- function(numSamples = 90) {
  out = matrix(nrow=4, ncol=numSamples)
  out[1, ] = rnorm(numSamples, 0, 1)
  out[2, ] = rnorm(numSamples, 0.2, 1)
  out[3, ] = rnorm(numSamples, 0, 1.2)
  out[4, ] = rnorm(numSamples, 0, 0.25)

  return(out)
}

generate.non.indep1 <- function(numSamples = 90) {
  M = 2
  out = matrix(nrow=2, ncol=numSamples)
  var = 0.5
  means = rnorm(n=3, 0, 10)
  w1 = runif(1)
  w2 = runif(1)
  for (i in 1:numSamples) {
    out[1, i] = normal.mixture.samples(c(means[1], means[2]), c(w1, 1 - w1), var, 1)
    out[2, i] = normal.mixture.samples(c(means[1], means[3]), c(w2, 1 - w2), var, 1)
  }
  return (out)
}

generate.non.indep2 <- function(numSamples = 90) {
  out = matrix(nrow=2, ncol=numSamples)
  for (i in 1:numSamples) {
    out[1, i] = rnorm(n=1, 0, 1)
    if (out[1, i] < 0)
      out[2, i] = rnorm(n=1, 1.5, 1)
    else
      out[2, i] = rnorm(n=1, -1.5, 1)
  }
  return (out)
}

get.distant.means <- function() {
  means = rnorm(n=3, 0, 10)
  if (abs(means[1] - means[2]) < 1 || abs(means[1] - means[3]) < 1)
    get.distant.means()
  else
    return(means)
}

generate.non.indep3 <- function(numSamples = 100) {
  out = matrix(nrow=2, ncol=numSamples)
  means = get.distant.means()
  w1 = runif(1, 0.25, 0.75)
  w2 = runif(1, 0.25, 0.75)
  var = 0.5
  out[1, ] = normal.mixture.samples(c(means[1], means[2]), c(w1, 1 - w1), var, numSamples)
  out[2, ] = normal.mixture.samples(c(means[1], means[3]), c(w2, 1 - w2), var, numSamples)
  return(out)
}
