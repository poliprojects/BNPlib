#include <iostream>
#include <string>
#include <vector>
#include <armadillo>
#include <sstream>


#include "../src/blocked_gibbs.hpp"
#include "../src/memory_collector_by_param.hpp"
#include "../src/memory_collector_by_iter.hpp"
#include "../src/normal_conjugate_hierarchy.hpp"
#include "../src/semi_hdp_blocked_gibbs.hpp"
#include "../src/marginal_semi_hdp.hpp"
#include "../src/utils.hpp"

std::vector<arma::vec> readData(std::string filename) {
  std::ifstream infile(filename);

  std::vector<std::vector<float>> out(3);

  int group;
  float grade;
  std::string line;
  char delim;
  // skip header
  std::getline(infile, line);
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    if (!(iss >> group >> delim >> grade) ) { break; }
    out[group - 1].push_back(grade);
  }

  std::vector<arma::vec> grades(3);
  for (int g=0; g < 3; g++) {
    grades[g] = arma::conv_to<arma::vec>::from(out[g]);
  }

  return grades;
};

void runexpmarginal(std::vector<arma::vec> data, std::string outdir) {
  int numIter = 100000;
  int thin = 10;
  int burnIn = 10000;

  std::string filename;

  MemoryCollectorByIter collector(numIter / thin);
  MarginalSemiHdp<NormalConjugateHierarchy> sampler(data);

  sampler.init();
  sampler.samplePseudoPrior(burnIn, burnIn);
  std::cout << "After pseudo prior" << std::endl;
  sampler.print();

  for (int i = 0; i < burnIn; i++) {
    sampler.sample();
  }

  for (int i = 0; i < numIter; i++) {
    sampler.sample();
    if (remainder(float(i), float(thin)) == 0) {
      collector.newIter();
      sampler.collect(&collector);
    }
  }

  std::cout << "\n\nFinal:" << std::endl;
  sampler.print();

  arma::mat predsMat = sampler.predictiveSamples(collector);
  filename = outdir + "marginal_preds.csv";
  predsMat.save(filename, arma::csv_ascii);
  filename = outdir + "marginal_c.csv";
  collector.getParamChain("c").save(filename, arma::csv_ascii);
}

int main() {
  std::string filename = "/home/mario/PhD/exchangeability/data/grades_chile_norm.csv";
  std::string outdir = "/home/mario/PhD/exchangeability/results/grades_chile/";

  std::vector<arma::vec> data = readData(filename);

  runexpmarginal(data, outdir);
}
