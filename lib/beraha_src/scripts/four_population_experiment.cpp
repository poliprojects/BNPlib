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
#include "../src/utils.hpp"


double rskewnorm(double alpha) {
  double u0 = stats::rnorm(0.0, 1.0);
  double v = stats::rnorm(0.0, 1.0);
  double d = alpha / sqrt(1 + alpha * alpha);
  double u1 = d * u0 + v * sqrt(1 - d * d);
  return  (u0 > 0) ? u1 : -u1;
}

arma::vec skewNormSample(int numSamples, double alpha) {
  arma::vec out(numSamples);
  for (int i=0; i < numSamples; i++)
    out(i) = rskewnorm(alpha);
  return out;
}

void runexpmarginal(std::vector<arma::vec> data, std::string outdir, int expnum) {
  int numIter = 100000;
  int thin = 10;
  int burnIn = 10000;
  // int numIter = 1000;
  // int thin = 10;
  // int burnIn = 1000;
  std::string filename;

  MemoryCollectorByIter collector(numIter / thin);
  MarginalSemiHdp<NormalConjugateHierarchy> sampler(data);

  sampler.init();
  sampler.samplePseudoPrior(burnIn, burnIn);
  std::cout << "After Pseudo Prior:";
  sampler.print();

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
  std::cout << "\n\nFinal: ";
  sampler.print();

  arma::mat predsMat = sampler.predictiveSamples(collector);
  filename = outdir + "marginal_exp" + std::to_string(expnum) + "_preds.csv";
  predsMat.save(filename, arma::csv_ascii);
  filename = outdir + "marginal_exp" + std::to_string(expnum) + "_c.csv";
  collector.getParamChain("c").save(filename, arma::csv_ascii);

  arma::vec xgrid = arma::linspace(-10, 10, 1000);
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
  std::string outdir = "/home/mario/PhD/exchangeability/results/four_pop/";

  std::vector<arma::vec> exp1Data(4);
  exp1Data[0] = arma::vec(200, arma::fill::randn);
  exp1Data[1] = arma::vec(200, arma::fill::randn);
  exp1Data[2] = arma::vec(200, arma::fill::randn);
  exp1Data[3] = skewNormSample(200, 1.0);
  experimentsData[0] = exp1Data;

  std::vector<arma::vec> exp2Data(4);
  exp2Data[0] = arma::vec(200, arma::fill::randn);
  exp2Data[1] = stats::rnorm<arma::vec>(200, 1.0, 0.0, 2.25);
  exp2Data[2] = stats::rnorm<arma::vec>(200, 1.0, 0.0, 0.25);
  exp2Data[3] = arma::vec(200, arma::fill::randn);
  experimentsData[1] = exp2Data;
  std::cout << "Data 2 ok" << std::endl;

  std::vector<arma::vec> exp3Data(4);
  exp3Data[0] = twoNormalMixture(200, 0.0, 1.0, 5.0, 1.0, 0.5);
  exp3Data[1] = twoNormalMixture(200, 0.0, 1.0, 5.0, 1.0, 0.5);
  exp3Data[2] = twoNormalMixture(200, 0.0, 1.0, -5.0, 1.0, 0.5);
  exp3Data[3] = twoNormalMixture(200, 5.0, 1.0, -5.0, 1.0, 0.5);
  experimentsData[2] = exp3Data;
  std::cout << "Data 3 ok" << std::endl;


  // runexpmarginal(experimentsData[2], outdir, 2 + 1);

  // #pragma omp parallel for
  for (uint i=0; i < experimentsData.size(); i++) {
    std::cout << "*** EXPERIMENT NUMBER : " << i << std::endl;
    runexpmarginal(experimentsData[i], outdir, i + 1);
    std::cout << "\n\n\n";
  }
};
