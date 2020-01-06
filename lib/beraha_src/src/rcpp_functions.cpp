#include "rcpp_functions.h"


Rcpp::List runGibbs(std::string model, const Rcpp::List& data_, int numIter,
                    int burnIn, int warmup, int thin, int numComponents,
                    int log_every, bool usePseudoPrior) {
  if (model == "normal_conjugate")
    return runNormalConjugateGibbs(
        data_, numIter, burnIn, warmup, thin, numComponents, log_every,
        usePseudoPrior);
  else if (model == "blocked_semi_hdp")
    return runSemiHdpBlockedGibbs(
        data_, numIter, burnIn, warmup, thin, numComponents, log_every,
        usePseudoPrior);
  else if (model == "marginal_semi_hdp")
    return runSemiHdpMarginal(
        data_, numIter, burnIn, warmup, thin, log_every, usePseudoPrior);
  else
    throw "Requested model is not implemented, supported are: "
          "'normal_conjugate', 'semi_hdp', 'marginal_semi_hdp'";
}

arma::mat predictiveSamples(std::string model, const Rcpp::List &chains) {
  arma::mat out;
  if (model == "blocked_semi_hdp") {
    MemoryCollectorByParam collector;
    collector.restoreFromR(chains);
    int numComponents = collector.get("clusterProbas", 0).n_cols;
    int numGroups = collector.get("clusterProbas", 0).n_rows;
    SemiHdpBlockedGibbs sampler(numGroups, numComponents);
    out = sampler.predictiveSamples(collector);
  } else if (model == "normal_conjugate") {
    MemoryCollectorByParam collector;
    collector.restoreFromR(chains);
    int numComponents = collector.get("clusterProbas", 0).n_cols;
    int numGroups = collector.get("clusterProbas", 0).n_rows;
    BlockedGibbs<NormalConjugateHierarchy> sampler(numGroups, numComponents);
    out = sampler.predictiveSamples(collector);
  } else if (model == "marginal_semi_hdp") {
    MemoryCollectorByIter collector;
    collector.restoreFromR(chains);
    int numGroups = collector.get("c", 0).n_elem;
    MarginalSemiHdp<NormalConjugateHierarchy> sampler(numGroups);
    out = sampler.predictiveSamples(collector);
  }
  return out;
}

Rcpp::List runNormalConjugateGibbs(
    const Rcpp::List& data_, int numIter, int burnIn, int warmup, int thin,
    int numComponents, int log_every, bool usePseudoPrior) {

  std::vector<arma::vec> data = Rcpp::as<std::vector<arma::vec>>(data_);
  BlockedGibbs<NormalConjugateHierarchy> sampler(data, numComponents);
  MemoryCollectorByParam collector(numIter / thin + 1);

  sampler.init();

  if (usePseudoPrior) {
    sampler.samplePseudoPrior(burnIn, warmup);
  } else {
    for (int i = 0; i < warmup; i++) {
      sampler.sample(true);
      if (remainder(double(i), double(log_every)) == 0) {
        Rcpp::Rcout << "Warm-up Iter: #" << i << " / " << warmup << std::endl;
        Rcpp::checkUserInterrupt();
      }
    }
  }

  for (int i = 0; i < burnIn; i++) {
    sampler.sample();
    if (remainder(double(i), double(log_every)) == 0) {
      Rcpp::Rcout << "Burn-in Iter: #" << i << " / " << burnIn << std::endl;
      Rcpp::checkUserInterrupt();
    }
  }

  for (int i = 0; i < numIter; i++) {
    sampler.sample();
    if (remainder(double(i + 1), double(thin)) == 0) {
      sampler.collect(&collector);
      Rcpp::checkUserInterrupt();
    }

    if (remainder(double(i), double(log_every)) == 0)
      Rcpp::Rcout << "Iter: #" << i << " / " << numIter << std::endl;
  }

  return collector.getChains();
}


Rcpp::List runSemiHdpBlockedGibbs(
    const Rcpp::List& data_, int numIter, int burnIn, int warmup, int thin,
    int numComponents, int log_every, bool usePseudoPrior) {

  std::vector<arma::vec> data = Rcpp::as<std::vector<arma::vec>>(data_);

  SemiHdpBlockedGibbs sampler(data, numComponents);
  MemoryCollectorByParam collector(numIter / thin + 1);

  sampler.init();

  if (usePseudoPrior) {
    sampler.samplePseudoPrior(burnIn, warmup);
    std::cout << "After Pseudo Prior" << std::endl;
    sampler.print();
  } else {
    for (int i = 0; i < warmup; i++) {
      sampler.sample(true);
      if (remainder(double(i), double(log_every)) == 0) {
        Rcpp::Rcout << "Warm-up Iter: #" << i << " / " << warmup << std::endl;
        Rcpp::checkUserInterrupt();
      }
    }
  }

  for (int i = 0; i < burnIn; i++) {
    sampler.sample();
    if (remainder(double(i), double(log_every)) == 0) {
      Rcpp::Rcout << "Burn-in Iter: #" << i << " / " << burnIn << std::endl;
      Rcpp::checkUserInterrupt();
    }
  }

  for (int i = 0; i < numIter; i++) {
    sampler.sample();
    if (remainder(double(i), double(thin)) == 0) {
      sampler.collect(&collector);
      Rcpp::checkUserInterrupt();
    }

    if (remainder(double(i), double(log_every)) == 0)
      Rcpp::Rcout << "Iter: #" << i << " / " << numIter << std::endl;
  }
  std::cout << "\n\nFinal" << std::endl;
  sampler.print();
  return collector.getChains();
}


Rcpp::List runSemiHdpMarginal(
    const Rcpp::List& data_, int numIter, int burnIn, int warmup, int thin,
    int log_every, bool usePseudoPrior) {
  Rcpp::Rcout << "NumIter: " << numIter << ", burnIn: "<< burnIn << ", warmup: " <<
    warmup << ", thin: " << thin << "log_every: " << log_every << ", usePseudoPrior: " <<
    usePseudoPrior << std::endl;

  std::vector<arma::vec> data = Rcpp::as<std::vector<arma::vec>>(data_);
  MarginalSemiHdp<NormalConjugateHierarchy> sampler(data);
  MemoryCollectorByIter collector(numIter / thin + 1);

  sampler.init();

  if (usePseudoPrior) {
    sampler.samplePseudoPrior(burnIn, warmup);
    Rcpp::Rcout << "Pseudoprior Done" << std::endl;
  } else {
    for (int i = 0; i < warmup; i++) {
      sampler.sample(true);
      if (remainder(double(i), double(log_every)) == 0) {
        Rcpp::Rcout << "Warm-up Iter: #" << i << " / " << warmup << std::endl;
        Rcpp::checkUserInterrupt();
      }
    }
  }

  for (int i = 0; i < burnIn; i++) {
    sampler.sample();
    if (remainder(double(i), double(log_every)) == 0) {
      Rcpp::Rcout << "Burn-in Iter: #" << i << " / " << burnIn << std::endl;
      Rcpp::checkUserInterrupt();
    }
  }

  for (int i = 0; i < numIter; i++) {
    sampler.sample();
    if (remainder(double(i), double(thin)) == 0) {
      collector.newIter();
      sampler.collect(&collector);
      Rcpp::checkUserInterrupt();
    }

    if (remainder(double(i), double(log_every)) == 0)
      Rcpp::Rcout << "Iter: #" << i << " / " << numIter << std::endl;
  }
  Rcpp::Rcout << "Finished MCMC Iterations" << std::endl;
  return collector.getChains();
}
