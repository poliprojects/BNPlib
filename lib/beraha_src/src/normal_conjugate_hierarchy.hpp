//
// Created by mario on 24/05/19.
//

#ifndef CPPMODEL_NORMAL_CONJUGATE_HIERARCHY_HPP
#define CPPMODEL_NORMAL_CONJUGATE_HIERARCHY_HPP

#include "options.hpp"
#include "utils.hpp"
#include "random_engine.hpp"


/*
 * This class represents a normal - normalgamma hierarchy
 */
class NormalConjugateHierarchy {
 public:
   class params {
    protected:
     double mu0, a, b, lambda; // state
     double mu0mean, mu0var, lambdaA, lambdaB;
     stats::rand_engine_t& engine = RandomEngine::Instance().get();
   public:
     params(
       double _a=5, double _b=5, double _mu0mean=0.0, double _mu0var=2.0,
       double _lambdaA=1.0, double _lambdaB=10.0);
     friend class NormalConjugateHierarchy;
     void update(arma::mat clusterVals);
     double getMu0() {return mu0;}
     double getA() {return a;}
     double getB() {return b;}
     double getLambda() {return lambda;}
   };

 protected:
  arma::vec param;
  std::shared_ptr<params> hyperparams;
  double postA, postB, postMean, postLambda;
  stats::rand_engine_t& engine = RandomEngine::Instance().get();
  int size = 2;

 public:

  ~NormalConjugateHierarchy() {}

  NormalConjugateHierarchy(params hyperparams) {
    this->hyperparams = std::make_shared<params>(hyperparams);
    param.resize(2, 1);
    samplePrior();
  }

  NormalConjugateHierarchy() {
    // HACK: should this thing be allowed?
    this->hyperparams = std::make_shared<params>(params());
    param.resize(2, 1);
    samplePrior();
  }

  NormalConjugateHierarchy(
      arma::vec param, params hyperparams):
       param(param) {
    this->hyperparams = std::make_shared<params>(hyperparams);
  }

  void samplePrior();

  void sample(const std::vector<double> &data);

  int getGTildeComponent() const {
    return -1;
  }

  const arma::vec &getParam() const {
    return param;
  }

  int getComponent() const {
    return 1;
  }

  double loglike(double x) {
    double sd = 1 / std::sqrt(param(1));
    return stats::dnorm(x, param(0), sd, true);
  }

  void reInitialize() {
    samplePrior();
  }

  double predict() {
    double sd = 1 / std::sqrt(param(1));
    return stats::rnorm(param(0), sd, engine);
  }

  void operator=(const NormalConjugateHierarchy& other) {
    this->param = other.param;
  }

  int getSize() {
    return size;
  }


  /*
   * Returns \int N(x | \mu, \sigma^2) dNIG(\mu, \sigma^2)
   */
  double marginalLogLike(double x) {
    arma::vec temp = normalGammaUpdate(
      arma::vec{x}, hyperparams->mu0, hyperparams->a, hyperparams->b, hyperparams->lambda);
    postMean = temp(0);
    postA = temp(1);
    postB = temp(2);
    postLambda = temp(3);
    double out = lgamma(postA) - lgamma(hyperparams->a);
    out += hyperparams->a * log(hyperparams->b) - postA * log(postB);
    out += 0.5 * (log(hyperparams->lambda) - log(postLambda));
    out -= M_PI;
    return out;
//    return stats::dt(x - mu0, 1, true);
  }

  static double marginalLogLike(const double x, const params& hypers) {
    arma::vec temp = normalGammaUpdate(
      arma::vec{x}, hypers.mu0, hypers.a, hypers.b, hypers.lambda);
    double postMean = temp(0);
    double postA = temp(1);
    double postB = temp(2);
    double postLambda = temp(3);
    double out = lgamma(postA) - lgamma(hypers.a);
    out += hypers.a * log(hypers.b) - postA * log(postB);
    out += 0.5 * (log(hypers.lambda) - log(postLambda));
    out -= M_PI;
    return out;
  }

};

#endif //CPPMODEL_NORMAL_CONJUGATE_HIERARCHY_HPP
