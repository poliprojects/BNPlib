//
// Created by mario on 24/05/19.
//

#include "normal_conjugate_hierarchy.hpp"


void NormalConjugateHierarchy::samplePrior() {
  param(1) = stats::rgamma(hyperparams->a, 1.0/hyperparams->b, engine);
  // TODO remove this line, fixing tau = 1 for debugging
  param(1) = 1.0;
  param(0) = stats::rnorm(
    hyperparams->mu0, std::sqrt(1.0 / (hyperparams->lambda * param(1))), engine);
}

void NormalConjugateHierarchy::sample(const std::vector<double> &data) {
  arma::vec datavec(data);
  //  std::cout << "Sampling for data: "; datavec.t().print();
  // normal inverse gamma update
  arma::vec temp = normalGammaUpdate(
    datavec, hyperparams->mu0, hyperparams->a, hyperparams->b, hyperparams->lambda);
  postMean = temp(0);
  postA = temp(1);
  postB = temp(2);
  postLambda = temp(3);
  // std::cout << "Post HyperParams: "; temp.t().print();
  double tau = stats::rgamma(postA, 1.0/postB, engine);
  // TODO: remove this line, fixing tau = 1 for debugging
//  tau = 1.0;
  double sigma = 1.0 / (std::sqrt(tau) * postLambda);
  double mu = stats::rnorm(postMean, sigma, engine);
  param(0) = mu;
  param(1) = tau;
  // std::cout << "mu: " << mu << ", tau: " << tau << std::endl;
}

NormalConjugateHierarchy::params::params(
  double _a, double _b, double _mu0mean, double _mu0var,
  double _lambdaA, double _lambdaB):
  a(_a), b(_b), mu0mean(_mu0mean), mu0var(_mu0var),
  lambdaA(_lambdaA), lambdaB(_lambdaB) {
    mu0 = mu0mean;
    lambda = _lambdaA / _lambdaB;
  }

void NormalConjugateHierarchy::params::update(arma::mat clusterVals) {
  double lambdaApost = lambdaA + 1.0 * clusterVals.n_rows / 2;
  arma::vec clusterMeans = clusterVals.col(0);
  arma::vec clusterVars = arma::pow(clusterVals.col(1), -1);
  double lambdaBPost = 0.5 * (
      lambdaB + arma::sum(arma::pow(clusterMeans - mu0mean, 2) / clusterVars));
  lambda = stats::rinvgamma(lambdaApost, 1.0 / lambdaBPost);
  double mu0precPost = \
    1.0 / mu0var + arma::sum(arma::pow(lambda * clusterVars, -1));
  double mu0meanPost = mu0mean / mu0var;
  for (int i=0; i<clusterMeans.n_elem; i++)
    mu0meanPost += clusterMeans(i) / (lambda * clusterVars(i));
  mu0meanPost /= mu0precPost;
  mu0 = stats::rnorm(mu0meanPost, 1 / std::sqrt(mu0precPost));
}
