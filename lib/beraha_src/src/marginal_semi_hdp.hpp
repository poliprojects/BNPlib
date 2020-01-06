//
// Created by mario on 27/06/19.
//

#ifndef CPPMODEL_MARGINAL_SEMI_HDP_HPP
#define CPPMODEL_MARGINAL_SEMI_HDP_HPP

#include "options.hpp"
#include "cluster.hpp"
#include "memory_collector_by_iter.hpp"
#include "random_engine.hpp"
#include "normal_conjugate_hierarchy.hpp"
#include "semi_hdp_hierarchy.hpp"
#include "utils.hpp"
#include "stirling_numbers.hpp"


template <typename Hierarchy>
class MarginalSemiHdp {
 protected:
  int numGroups;
  arma::uvec samplesPerGroup;
  std::vector<arma::vec> data;

  bool usePseudoPrior = false;
  MemoryCollectorByIter pseudoPriorCollector;
  int numPseudoPriorIters;

  arma::uvec c;
  std::vector<std::vector<int>> datum2cluster;
  std::vector<std::vector<bool>> datum2hierarchy;
  std::vector< std::vector<Cluster<Hierarchy>> > baseClusters;
  std::vector<Cluster<Hierarchy>> hdpClusters;
  arma::vec betas;
  arma::vec numTables;

  typename Hierarchy::params G0Params;
  typename Hierarchy::params G00Params;

  std::vector< std::vector<Cluster<Hierarchy>> > pseudoBaseClusters;
  std::vector<Cluster<Hierarchy>> pseudoHdpClusters;

  arma::uvec nDataToHdpFromRestaurant;
  int nDataToHdp;

  double semiHdpWeight;
  double alpha; // concentration of F_i DPs
  double gamma; // concentration of Gtilde DP
  double alpha_0; // parameter of rho's Polya-Urn
  double w_a, w_b; // parameters for w's prior [w ~ Beta(w_a, w_b)]
  arma::vec omegas; // parameters for c's prior [c_i ~ Cat(omegas)]
  arma::vec etas; // parameters for omegas prior [omegas ~ Dirichlet(etas)]

  stats::rand_engine_t& engine = RandomEngine::Instance().get();

 public:
  ~MarginalSemiHdp() = default;
  MarginalSemiHdp() {}

  MarginalSemiHdp(std::vector<arma::vec> & _data);

  MarginalSemiHdp(int numGroups): numGroups(numGroups) {
    baseClusters.resize(numGroups);
  }

  void init();

  void sample(bool warmpup=false) {
    sampleAssignments();
    sampleHieararchy();
    sampleClusters();
    if (! warmpup)
      sampleGroupAssignments();
    relabel();
    sampleSemiHdpWeight();
    sampleHypers();
  }

  void samplePseudoPrior(int burin, int numsteps);

  void sampleClusters();

  void sampleAssignments();

  void sampleAssignmentsForGroup(int group, bool flush=false);

  void sampleHieararchy();

  void relabel();

  void sampleGroupAssignments();

  void sampleSemiHdpWeight();

  void sampleHdpLatentNtables();

  void sampleHdpLatentBetas();

  void sampleHypers();

  double logLikeForComponent(double x, int component);

  void collect(MemoryCollectorByIter* collector);

  arma::mat predictiveSamples(const MemoryCollectorByIter& collector);

  double _probaForHdpClus(int h, double x, int component) {
    int n0rh = hdpClusters[h].size(component);
    double beta_h = 0;
    if (h < betas.n_elem)
      beta_h = betas(h);
    else {
      beta_h = 0.0;
      std::cout << "PD ";
      assert(-1 > 0);
    }

    // double beta_h = (h < betas.n_elem) ? betas(h) : 1.0;
    double proba = 1.0 * n0rh + alpha * beta_h;
    double loglike = hdpClusters[h].getParam().loglike(x);
    proba *= exp(loglike);

    return proba;
  }

  double _probaForNewHdpClus(double x, int component) {
//    double den = (nDataToHdpFromRestaurant[component] -1 + alpha);
    double logLike = Hierarchy::marginalLogLike(x, G00Params);
    double out = alpha * betas(betas.n_elem - 1) * exp(logLike);

    return out;
  }

  int _assignDatumToRestaurant(double datum, int group, int restaurant, bool addNew=true);

  int _assignDatumToHdp(double datum, int group, int restaurant, bool addNew=true);

  void _restoreState(const MemoryCollectorByIter& collector, int iternum);

  void print();

  std::vector<arma::mat> estimateDensities(
      const arma::vec& grid, const MemoryCollectorByIter& collector) {
    std::vector<arma::vec> grids(numGroups);
    for (int i=0; i < numGroups; i++)
      grids[i] = grid;

    return estimateDensities(grids, collector);
  }

  std::vector<arma::mat> estimateDensities(
      const std::vector<arma::vec>& grids,
      const MemoryCollectorByIter& collector);

  const arma::uvec &getC() const;
  const std::vector< std::vector<Cluster<Hierarchy>> > getBaseClusters() const {
    return baseClusters;
  }
};

#include "marginal_semi_hdp_imp.hpp"

#endif //CPPMODEL_MARGINAL_SEMI_HDP_HPP
