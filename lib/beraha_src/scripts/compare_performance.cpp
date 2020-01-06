#include <chrono>
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
#include "../src/utils.hpp"

int numIter = 1000;
int thin = 10;
int burnIn = 1000;

void runexpblocked(std::vector<arma::vec> data, int numComponents) {
  BlockedGibbs<NormalConjugateHierarchy> sampler(data, numComponents);

  sampler.init();
  sampler.samplePseudoPrior(burnIn, burnIn);

  for (int i = 0; i < burnIn; i++) {
    sampler.sample();
  }

  for (int i = 0; i < numIter; i++) {
    sampler.sample();
  }
}


void runexpmarginal(std::vector<arma::vec> data) {
  MarginalSemiHdp<NormalConjugateHierarchy> sampler(data);

  sampler.init();
  sampler.samplePseudoPrior(burnIn, burnIn);

  for (int i = 0; i < burnIn; i++) {
    sampler.sample();
  }

  for (int i = 0; i < numIter; i++) {
    sampler.sample();
  }
}

int main() {
  std::string outfile = "/home/mario/PhD/exchangeability/results/two_pop/performance_comparison.csv";

  std::vector<int> numComponents{10, 20, 40, 80};
  std::vector<int> numSamples{100, 200, 400, 1000, 2000};
  arma::mat out(numSamples.size(), numComponents.size() + 1);

  arma::vec means(10);
  arma::vec sds(10, arma::fill::eye);
  for (int i=0; i < 10; i++) {
    means(i) = stats::rnorm(0.0, 10.0);
  }

  std::vector<arma::vec> data(2);
  for (int row=0; row < numSamples.size(); row++) {
    data[0] = uniformNormalMixture(numSamples[row], means, sds);
    data[1] = uniformNormalMixture(numSamples[row], means, sds);
    for (int col=0; col < numComponents.size(); col++) {
      std::cout << "Num Samples: " << numSamples[row] <<
                   ", Num Components: " << numComponents[col] << std::endl;
      auto start = std::chrono::steady_clock::now();
      runexpblocked(data, numComponents[col]);
      auto end = std::chrono::steady_clock::now();
      out(row, col) = std::chrono::duration<double, std::milli> (end - start).count();
    }
    auto startMarg = std::chrono::steady_clock::now();
    std::cout << "Starting Marginal Sampler" << std::endl;
    runexpmarginal(data);
    std::cout << "Ended Marginal Sampler" << std::endl;
    auto endMarg = std::chrono::steady_clock::now();
    out(row, numComponents.size()) = \
      std::chrono::duration<double, std::milli> (endMarg - startMarg).count();
  }
  out.save(outfile, arma::csv_ascii);
}
