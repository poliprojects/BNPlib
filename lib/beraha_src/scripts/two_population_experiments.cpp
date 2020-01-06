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


void runexpblocked(std::vector<arma::vec> data, std::string outdir, int expnum) {
  int numIter = 100000;
  int thin = 10;
  int burnIn = 10000;
  // int numIter = 1000;
  // int thin = 10;
  // int burnIn = 1000;
  int numComponents = 20;
  std::string filename;

  MemoryCollectorByParam collector(numIter / thin);
  SemiHdpBlockedGibbs sampler(data, numComponents);

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

  arma::mat predsMat = sampler.predictiveSamples(collector);
  filename = outdir + "blocked_exp" + std::to_string(expnum) + "_preds.csv";
  predsMat.save(filename, arma::csv_ascii);
  filename = outdir + "blocked_exp" + std::to_string(expnum) + "_c.csv";
  collector.getParamChain("c").save(filename, arma::csv_ascii);
  filename = outdir + "blocked_exp" + std::to_string(expnum) + "_chains.csv";
  collector.to_csv(filename);
}


void runexpmarginal(std::vector<arma::vec> data, std::string outdir, int expnum) {
  int numIter = 10000;
  int thin = 10;
  int burnIn = 100000;
  // int numIter = 1000;
  // int thin = 10;
  // int burnIn = 1000;
  std::string filename;

  MemoryCollectorByIter collector(numIter / thin);
  MarginalSemiHdp<NormalConjugateHierarchy> sampler(data);

  sampler.init();
  sampler.samplePseudoPrior(burnIn, burnIn);

  for (int i = 0; i < burnIn; i++) {
    sampler.sample();
  }

  for (int i = 0; i < numIter; i++) {
    sampler.sample();
    if (remainder(double(i), double(thin)) == 0) {
      collector.newIter();
      sampler.collect(&collector);
    }
  }

  arma::mat predsMat = sampler.predictiveSamples(collector);
  filename = outdir + "marginal_exp" + std::to_string(expnum) + "_preds.csv";
  predsMat.save(filename, arma::csv_ascii);
  filename = outdir + "marginal_exp" + std::to_string(expnum) + "_c.csv";
  collector.getParamChain("c").save(filename, arma::csv_ascii);
  std::cout << "Finished experiment: " << expnum << std::endl;

  arma::vec xgrid = arma::linspace(-10, 15, 1000);
  std::vector<arma::mat> densities = sampler.estimateDensities(xgrid, collector);
  for (int i=0; i < data.size(); i++) {
    filename = outdir + "marginal_exp" + std::to_string(expnum) + \
              + "_group" + std::to_string(i+1) + "_estimated_dens.csv";
    densities[i].save(filename, arma::csv_ascii);
  }
}


int main() {

  std::vector<std::vector<arma::vec>> experimentsData(3);
  // std::vector<std::vector<arma::vec>> experimentsData(1);
  std::string outdir = "/home/mario/PhD/exchangeability/results/two_pop/";

  std::vector<arma::vec> exp1Data(2);
  // exp1Data[0] = twoNormalMixture(100, 0.0, 1.0, 5.0, 1.0, 0.5);
  // exp1Data[1] = twoNormalMixture(100, 0.0, 1.0, 5.0, 1.0, 0.5);
  exp1Data[0] = arma::vec(100, arma::fill::randn);
  exp1Data[0] = arma::join_cols(exp1Data[0], arma::vec(100, arma::fill::randn) + 5);
  exp1Data[0] = arma::shuffle(exp1Data[0]);
  exp1Data[1] = arma::vec(100, arma::fill::randn);
  exp1Data[1] = arma::join_cols(exp1Data[1], arma::vec(100, arma::fill::randn) + 5);
  exp1Data[1] = arma::shuffle(exp1Data[1]);
  experimentsData[0] = exp1Data;

  std::vector<arma::vec> exp2Data(2);
  exp2Data[0] = twoNormalMixture(100, 5.0, 0.6, 10.0, 0.6, 0.9);
  exp2Data[1] = twoNormalMixture(100, 5.0, 0.6, 0.0, 0.6, 0.1);
  experimentsData[1] = exp2Data;

  std::vector<arma::vec> exp3Data(2);
  exp3Data[0] = twoNormalMixture(100, 5.0, 1.0, 0.0, 1.0, 0.8);
  exp3Data[1] = twoNormalMixture(100, 5.0, 1.0, 0.0, 1.0, 0.2);
  experimentsData[2] = exp3Data;

  // std::vector<arma::vec> exp4Data(2);
  // exp4Data[0] = twoNormalMixture(100, 5.0, 0.6, 10.0, 0.6, 0.5);
  // exp4Data[1] = twoNormalMixture(100, 5.0, 0.6, 0.0, 0.6, 0.5);
  // experimentsData[3] = exp4Data;

  std::cout << "Starting" << std::endl;

  // #pragma omp parallel for
  for (uint i=0; i < experimentsData.size(); i++) {
    std::cout << "i: " << i << std::endl;
    // runexpblocked(experimentsData[i], outdir, i + 1);
    runexpmarginal(experimentsData[i], outdir, i + 1);
  }
};
