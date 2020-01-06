#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <armadillo>


#include "../src/blocked_gibbs.hpp"
#include "../src/memory_collector_by_param.hpp"
#include "../src/memory_collector_by_iter.hpp"
#include "../src/normal_conjugate_hierarchy.hpp"
#include "../src/semi_hdp_blocked_gibbs.hpp"
#include "../src/marginal_semi_hdp.hpp"
#include "../src/metropolised_marginal_semi_hdp.hpp"
#include "../src/utils.hpp"




int main() {
  int N = 100;

  std::vector<arma::vec> data(N);
  // std::vector<std::vector<arma::vec>> experimentsData(1);
  std::string outdir = "/home/mario/PhD/exchangeability/results/many_pop/";

  for (int i = 0; i < N / 2; i++) {
    data[i] = twoNormalMixture(200, 0.0, 1.0, 5.0, 1.0, 0.5);
    data[N/2 + i] = twoNormalMixture(200, 0.0, 1.0, -5.0, 1.0, 0.5);
  }


  int numIter = 10000;
  int thin = 10;
  int burnIn = 10000;

  std::string filename;

  MemoryCollectorByIter collector(numIter / thin);
  MetropolisedMarginalSemiHdp<NormalConjugateHierarchy> sampler(data);

  sampler.init();
  sampler.samplePseudoPrior(burnIn, burnIn);
  std::cout << "Finished Pseudo Prior" << std::endl;

  for (int i = 0; i < burnIn; i++) {
    sampler.sample();

    if (remainder(double(i), double(200)) == 0) {
      std::cout << "Burn-In Iter: " << i << " / " << burnIn << std::endl;
    }
  }

  for (int i = 0; i < numIter; i++) {
    sampler.sample();
    if (remainder(double(i), double(thin)) == 0) {
      collector.newIter();
      sampler.collect(&collector);
    }

    if (remainder(double(i), double(100)) == 0) {
      std::cout << "Iter: " << i << " / " << numIter << std::endl;
    }
  }

  arma::mat predsMat = sampler.predictiveSamples(collector);
  filename = outdir + "marginal_preds.csv";
  predsMat.save(filename, arma::csv_ascii);
  filename = outdir + "marginal_c.csv";
  collector.getParamChain("c").save(filename, arma::csv_ascii);
};
