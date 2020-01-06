#include <iostream>
#include <string>
#include <vector>
#include <armadillo>


#include "../src/blocked_gibbs.hpp"
#include "../src/memory_collector_by_param.hpp"
#include "../src/normal_conjugate_hierarchy.hpp"
#include "../src/semi_hdp_blocked_gibbs.hpp"
#include "../src/utils.hpp"


void runexp(std::vector<arma::vec> data, std::string outdir, int expnum) {
  int numIter = 100000;
  int thin = 10;
  int burnIn = 10000;
  // int numIter = 1000;
  // int thin = 10;
  // int burnIn = 1000;
  int numComponents = 20;
  std::string filename;

  MemoryCollectorByParam collector(numIter / thin);
  BlockedGibbs<NormalConjugateHierarchy> sampler(data, numComponents);
  sampler.init();
  sampler.samplePseudoPrior(burnIn, burnIn);

  for (int i = 0; i < burnIn; i++) {
    sampler.sample();
  }

  for (int i = 0; i < numIter; i++) {
    sampler.sample();
    if (remainder(double(i), double(thin)) == 0)
      sampler.collect(&collector);
  }

  filename = outdir + "exp" + std::to_string(expnum) + "_c_blocked.csv";
  collector.getParamChain("c").save(filename, arma::csv_ascii);
}

int main() {

  std::string outdir = "/home/mario/PhD/exchangeability/results/blocked/non_indep/";

  std::cout << "Starting" << std::endl;

  #pragma omp parallel for
  for (uint i=0; i < 50; i++) {
    double mean1, mean2, mean3, sd1, sd2, sd3, w;
    std::vector<arma::vec> data(2);
    mean1 = stats::rnorm(0.0, 10.0);
    mean2 = stats::rnorm(0.0, 10.0);
    mean3 = stats::rnorm(0.0, 10.0);
    sd1 = stats::rgamma(2.0, 0.5);
    sd2 = stats::rgamma(2.0, 0.5);
    sd3 = stats::rgamma(2.0, 0.5);
    w = stats::rbeta(1.0, 1.0);
    data[0] = twoNormalMixture(100, mean1, sd1, mean2, sd2, w);
    data[1] = twoNormalMixture(100, mean1, sd1, mean3, sd3, w);
    std::cout << "i: " << i << std::endl;
    runexp(data, outdir, i);
  }
};
