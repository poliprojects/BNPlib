//
// Created by mario on 13/05/19.
//

#ifndef CPPMODEL_SEMI_HDP_HIERARCHY_HPP
#define CPPMODEL_SEMI_HDP_HIERARCHY_HPP

#include "options.hpp"
#include "random_engine.hpp"
#include "utils.hpp"
#include "metropolis_hastings.hpp"


/*
 * Sampler for the parameter theta = (mu, sigma^2)
 * We perform a change of variables: sigma^2 -> eta = log(sigma^2) not to deal
 * with constrained distributions.
 * The state is (mu, eta), so the exp transformation must be done afterwards.
 */
class ThetaSampler: public MetropolisHastings<arma::vec> {
 protected:
  arma::vec data;
  double mu0;
  double lambda;
  double a;
  double b;

  std::shared_ptr<arma::mat> gTildeAtoms;
  std::shared_ptr<arma::vec> gTidleProbs;
  std::shared_ptr<double> gTildeWeight;
  int gTildeComp;

 public:
  ~ThetaSampler() {gTildeComp = -2;}

  ThetaSampler(): MetropolisHastings(arma::vec{0.0, 0.0}) {gTildeComp = -2;}

  ThetaSampler(
      std::shared_ptr<arma::mat> const& _gTidleAtoms,
      std::shared_ptr<arma::vec> const& _gTidleProbs,
      std::shared_ptr<double> _gTildeWeight):
        MetropolisHastings(arma::vec{0.0, 1.0}),
        gTildeAtoms(_gTidleAtoms), gTidleProbs(_gTidleProbs),
        gTildeWeight(_gTildeWeight) {
    gTildeComp = -2;
  }

  /* Prior = wG0 + (1-w) Gtilde
   * G0: N(mu | 0, sigma^2 * lambda) x InvGamma(sigma^2 | a, b)
   * Change of variables: eta = log(sigma^2)
   * -> G0: N(mu | 0, exp(eta) * lambda) x InvGamma(exp(eta) | a, b) x exp(eta)
   *
   * state  is (mu, eta)
   */
  double logprior(arma::vec state) {
    double g0contrib = stats::dnorm(
        state(0), mu0, std::sqrt(exp(state(1)) * (lambda)), true);
    g0contrib += stats::dinvgamma(exp(state(1)), a, b, true);
    g0contrib += state(1);
    g0contrib += log(*gTildeWeight);

    double gTildeContrib = EPS;
    arma::rowvec constrainedState = currState.t();
    constrainedState(1) = exp(constrainedState(1));
    for (int i = 0; i < gTildeAtoms->n_rows; i++) {
      if (arma::all(constrainedState - gTildeAtoms->row(i) < EPS))
        gTildeContrib += log(gTidleProbs->at(i));
    }
    gTildeContrib += log(1- *gTildeWeight);

    return logSumExp(arma::vec{g0contrib, gTildeContrib});
  }

  double loglike(arma::vec state) {
    // TODO ACK: should probably call SemiHdpHierarchy's method...
    double out=0.0;
    for (int i = 0; i < data.n_elem; i++)
      out += stats::dnorm(data[i], state(0), std::sqrt(exp(state(1))), true);
    return out;
  }

  double logpost(arma::vec state) {
    return logprior(state) + loglike(state);
  }

  double logratio(arma::vec newState) {
    return logpost(newState) - logpost(currState);
  }

  void setData(const arma::vec& data_) {
    data = data_;
  }

  void setHyperParams(double _mu0, double _lambda, double _a, double _b) {
    mu0 = _mu0;
    lambda = _lambda;
    a = _a;
    b= _b;
  }

  /*
   * Very silly proposal: with probability 0.5 we propose a draw from a
   * bivariate normal for (mu, eta), else we draw uniformly from gTilde's atoms
   */
  arma::vec propose() {
    double r = stats::runif(0, 1, engine);
    arma::vec out(2);
    if (r < *gTildeWeight) {
      out(0) = stats::rnorm(currState(0), 0.1, engine);
      out(1) = stats::rnorm(currState(1), 0.1, engine);
      gTildeComp = -1;
    } else {
      arma::rowvec constrainedState = currState.t();
      constrainedState(1) = exp(constrainedState(1));
      arma::vec probas(gTildeAtoms->n_rows);
      for (int i = 0; i < gTildeAtoms->n_rows; i++) {
        probas(i) = 1.0 / log(
            arma::norm(gTildeAtoms->row(i) - constrainedState));
      }
      probas /= arma::sum(probas);
      int row = stats::rdiscrete(probas, engine);
      assert(row != gTildeAtoms->n_rows);
      out = gTildeAtoms->row(row).t();
      out(1) = log(out(1));
//      std::cout << "Curr state: "; currState.t().print();
//      std::cout << "Proposing from gTilde: "; out.t().print();
      gTildeComp = row;
    }
    return out;
  }

  int getGTildeComp() const {
    return gTildeComp;
  }

};

/*
 * Hierachy of the Semi HDP
 * Y_ij | c_i, s_ij, means, vars ~ Normal(means[c_i, s_ij], vars[c_i, s_ij])
 * means, vars ~ w G_0 + (1 - w) Gtilde
 * Gtilde ~ DP(gamma, P_00)
 *
 * G_0 = NormalInvGamma
 * P_00 = NormalInvGamma
 *
 * We augment the state with the variable comp ~ Bern(w), so that if
 * comp = 1, the update is a Normal-NormalInvGamma, if comp = 0,
 * the update selects among the atoms of Gtilde
 *
 * This class can sample a mean and variance pair from the prior and
 * perform a metropolis step for the posterior.
 */

class SemiHdpHierarchy {
 protected:
  // (mean, var)
  arma::vec param;

  std::shared_ptr<double> w;
  int comp;

  // Matrix (numcomponents x 2). It is the realization of Gtilde
  std::shared_ptr<arma::mat> gTildeAtoms;
  std::shared_ptr<arma::vec> gTildeProbs;
  int gTildeComponent;

//  double mu0, a, b, lambda;
  double mu0 = 0.0;
  double a = 2.0;
  double b = 2.0;
  double lambda = 1.0;
  double postA, postB, postMean, postLambda;
  int numComponents;
  int size = 2;

  ThetaSampler thetaSampler;
  stats::rand_engine_t& engine = RandomEngine::Instance().get();

 public:

  ~SemiHdpHierarchy() {}

  SemiHdpHierarchy() {}

  SemiHdpHierarchy(arma::vec param): param(param) {}

  SemiHdpHierarchy(
      std::shared_ptr<double> const& _w,
      std::shared_ptr<arma::mat> const& _gTildeAtoms,
      std::shared_ptr<arma::vec> const& _gTidleProbs):
        w(_w), gTildeAtoms(_gTildeAtoms), gTildeProbs(_gTidleProbs),
        thetaSampler(_gTildeAtoms, _gTidleProbs, _w) {
    param.resize(2, 1);

    numComponents = gTildeProbs->n_elem;

    // FIXME Fixing Hyperparams
    mu0 = 0.0;
    a = 1.0;
    b = 1.0;
    lambda = 1.0;

    comp = stats::rbern(*w);
    samplePrior();

    gTildeComponent = -2;
    thetaSampler.setHyperParams(mu0, lambda, a, b);
  }

  void samplePrior() {
    // TODO HERE
    param(1) = 1 / stats::rgamma(a, 1 / b, engine);
    param(0) = stats::rnorm(mu0, param(1) / lambda, engine);
    gTildeComponent = -1;
//    if (stats::rbern(*w, engine) == 1) {
//      param(1) = stats::rinvgamma(a, b, engine);
//      param(0) = stats::rnorm(mu0, param(1) / lambda, engine);
//    } else {
//      int rowNum = stats::rdiscrete(*gTildeProbs, engine);
//      param(0) = gTildeAtoms->row(rowNum)(0);
//      param(1) = gTildeAtoms->row(rowNum)(1);
//    }
  }

  void sample(const std::vector<double> &data) {
    arma::vec datavec(data);
    arma::vec state = param;
    state(1) = log(state(1));
    thetaSampler.setCurrState(state);
    thetaSampler.setData(datavec);
    metropolisHastingsUpdate();
  }

  void metropolisHastingsUpdate() {
    thetaSampler.step();
    gTildeComponent = thetaSampler.getGTildeComp();
    if (gTildeComponent > -1)
      comp = 0;
    else
      comp = 1;
    param = thetaSampler.getState();
    param(1) = exp(param(1));
  }

  void computeNormInvGammaUpdate(const arma::vec& datavec) {
    // normal inverse gamma update
    arma::vec temp = normalGammaUpdate(datavec, mu0, a, b, lambda);
    postMean = temp(0);
    postA = temp(1);
    postB = temp(2);
    postLambda = temp(3);
  }

  void sampleParam(const arma::vec& datavec) {
    int n = datavec.size();
    if (comp == 1) {
      double sigma = 1 / stats::rgamma(postA, 1 / postB, engine);
      double mu = stats::rnorm(postMean, std::sqrt(sigma / postLambda), engine);
      param(0) = mu;
      param(1) = sigma;
      gTildeComponent = -1;
    } else {
      // choose the atom from Gtildes ones
      arma::vec probas(gTildeAtoms->n_rows);
      for (int h = 0; h < numComponents; h++) {
        double mean = gTildeAtoms->row(h)[0];
        double sd = std::sqrt(gTildeAtoms->row(h)[1]);
        arma::vec loglikes(datavec.n_elem);
        for (int i = 0; i < datavec.n_elem; i++)
          loglikes(i) = stats::dnorm(datavec(i), mean, sd, true);

        probas(h) = exp(log(gTildeProbs->at(h)) + arma::sum(loglikes));
      }

      probas /= arma::sum(probas);
      gTildeComponent = stats::rdiscrete(probas, engine);
      param(0) = gTildeAtoms->row(gTildeComponent)(0);
      param(1) = gTildeAtoms->row(gTildeComponent)(1);
    }
  }

  void sampleComp(const arma::vec& datavec) {
    arma::vec logprobas(2);
    double n = datavec.n_elem;
    logprobas(0) = log(*w) + \
      lgamma(postA) - lgamma(a) + a * log(b) - postA * log(postB) + \
      0.5 * (log(postLambda) - log(lambda)) - n/2 * log(2 * M_PI);
//    logprobas(0) = log(*w) + \
//      stats::dnorm(param(0), mu0, param(1) / lambda, true) +
//      stats::dinvgamma(param(1), a, b, true);


    logprobas(1) = log(1 - *w);
    for (int i = 0; i < n; i++) {
      arma::vec curr(numComponents);
      for (int h = 0; h < numComponents; h++) {
        arma::vec param = gTildeAtoms->row(h).t();
        curr(h) = log(gTildeProbs->at(h)) + \
                  stats::dnorm(datavec(i), param(0),
                      std::sqrt(param(1)), true);
      }
      logprobas(1) += logSumExp(curr);
    }
    arma::vec probas = arma::exp(logprobas);
    probas /= arma::sum(probas);
    comp = stats::rbern(probas(0), engine);
  }

  void reInitialize() {
    samplePrior();
    gTildeComponent = -1;
  }

  double predict() {
    return stats::rnorm(param(0), std::sqrt(param(1)), engine);
  }

  const arma::vec &getParam() const {
    return param;
  }

  double loglike(double x) {
    return stats::dnorm(x, param(0), std::sqrt(param(1)), true);
  }

  static double loglike(double x, const arma::vec& param) {
    return stats::dnorm(x, param(0), std::sqrt(param(1)), true);
  }

  int getComponent() const {
    return comp;
  }

  void setComp(int comp) {
    SemiHdpHierarchy::comp = comp;
  }

  int getGTildeComponent() const {
    return gTildeComponent;
  }

  void operator=(const SemiHdpHierarchy &other) {
    this->param = other.param;
    thetaSampler.setCurrState(this->param);
  }

  int getSize() {
    return size;
  }
};

#endif //CPPMODEL_SEMI_HDP_HIERARCHY_HPP
