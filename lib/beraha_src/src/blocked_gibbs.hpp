//
// Created by mario on 30/04/19.
//

#ifndef CPPMODEL_NORMAL_UNIV_HPP
#define CPPMODEL_NORMAL_UNIV_HPP

#include "options.hpp"
#include "cluster.hpp"
#include "memory_collector_by_param.hpp"
#include "random_engine.hpp"
#include "normal_conjugate_hierarchy.hpp"
#include "semi_hdp_hierarchy.hpp"
#include "utils.hpp"


/*
 * Perform Gibbs sampling for univariate normal DPMs
 */
template <typename  Hierarchy>
class BlockedGibbs {
 protected:
  int numGroups;
  int numComponents;
  arma::uvec samplesPerGroup;
  std::vector<arma::vec> data;

  bool usePseudoPrior = false;
  MemoryCollectorByParam pseudoPriorCollector;
  int numPseudoPriorIters;

  // theta_{ij} = theta*{lm} <--> c[i] = l, datum2cluster[i][j] = m
  arma::uvec c;
  std::vector<std::vector<int>> datum2cluster;
  // TODO here
  std::vector< std::vector<Cluster<Hierarchy>> > clusters;
  arma::mat clusterProbas;

  double alpha; // concentration of F_i DPs
  double gamma; // concentration of Gtilde DP

  double alpha_0; // parameter of rho's Polya-Urn

  stats::rand_engine_t& engine = RandomEngine::Instance().get();

 public:
  ~BlockedGibbs() = default;

  BlockedGibbs(int numGroups, int H);

  BlockedGibbs(std::vector<arma::vec> & _data, int H);

  void init();

  void sample(bool warmpup=false) {
    sampleClusters();
    sampleClusterProbs();
    sampleAssignments();
    if (! warmpup)
      sampleGroupAssignments();
  }

  void samplePseudoPrior(int burnin, int numSteps);

  arma::uvec samplePolyaUrn(double concentration, int size);

  void sampleClusters();

  void sampleAssignments();

  void sampleAssignmentsForGroup(int group, bool flush=false);

  void sampleClusterProbs();

  arma::vec _updateClusterProbas(arma::vec sizes, double alpha);

  void sampleGroupAssignments();

  double logLikeForComponent(double x, int component);

  void collect(MemoryCollectorByParam* collector);

  arma::mat predictiveSamples(const MemoryCollectorByParam& collector);

  void print();

  int getNumGroups() const;
  const arma::uvec &getSamplesPerGroup() const;
  const arma::uvec &getC() const;
  const arma::mat &getClusterProbas() const;
};

#include "blocked_gibbs_imp.hpp"

#endif //CPPMODEL_NORMAL_UNIV_HPP
