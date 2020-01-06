//
// Created by mario on 16/05/19.
//

#include "utils.hpp"

double logSumExp(arma::vec logprobas) {
  double maxVal = arma::max(logprobas);
  double out = 0.0;

  for (int i=0; i < logprobas.n_elem; i++) {
    out += exp(logprobas(i) - maxVal);
  }

  return log(out) + maxVal;

}

arma::vec normalGammaUpdate(
    arma::vec data, double priorMean, double priorA, double priorB,
    double priorLambda) {
  double postMean, postA, postB, postLambda;
  int n = data.size();
  if (n == 0) {
    return arma::vec{priorMean, priorA, priorB, priorLambda};
  }
  double ybar = arma::mean(data);
  postMean = (priorLambda * priorMean + n * ybar) / (priorLambda + n);
  postA = 1.0 * priorA + 1.0 * n / 2;

  // arma::var(x, 1) divides by n, not n-1
  double ss = n * arma::var(data, 1);

  postB = (
      priorB + 0.5 * ss +
      0.5 * priorLambda / (n + priorLambda) * n * std::pow((ybar - priorMean), 2));

  postLambda = priorLambda + n;

  return arma::vec{postMean, postA, postB, postLambda};
}


arma::vec twoNormalMixture(int numSamples, double mean1, double sd1,
                           double mean2, double sd2, double w) {
  arma::vec out(numSamples);
  double mean, sd;
  for (int i=0; i < numSamples; i++) {
    int comp1 = stats::rbern(w);
    if (comp1 == 1) {
      mean = mean1;
      sd = sd1;
    } else {
      mean = mean2;
      sd = sd2;
    }
    out(i) = stats::rnorm(mean, sd);
  }
  return out;
}

arma::vec uniformNormalMixture(int numSamples, arma::vec means, arma::vec sds) {
  arma::vec out(numSamples);
  uint numComp = means.n_elem;
  for (int i=0; i < numSamples; i++) {
    int comp = stats::rdiscreteunif(0, numComp);
    out(i) = stats::rnorm(means[comp], sds[comp]);
  }
  return out;
}
