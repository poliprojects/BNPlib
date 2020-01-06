//
// Created by mario on 27/06/19.
//

#ifndef CPPMODEL_METROPOLISED_MARGINAL_SEMI_HDP_HPP
#define CPPMODEL_METROPOLISED_MARGINAL_SEMI_HDP_HPP

#include "options.hpp"
#include "cluster.hpp"
#include "memory_collector_by_iter.hpp"
#include "random_engine.hpp"
#include "normal_conjugate_hierarchy.hpp"
#include "semi_hdp_hierarchy.hpp"
#include "utils.hpp"
#include "stirling_numbers.hpp"
#include "marginal_semi_hdp.hpp"


template <typename Hierarchy>
class MetropolisedMarginalSemiHdp: public MarginalSemiHdp<Hierarchy> {
 protected:
   arma::mat distances;

 public:

  ~MetropolisedMarginalSemiHdp() {}
  MetropolisedMarginalSemiHdp(std::vector<arma::vec> & _data):
      MarginalSemiHdp<Hierarchy>(_data) {
    distances.resize(numGroups, numGroups);
    distances.zeros();
  }

  /*
  * Computes L2 distance between two mixture of gaussians
  */
  double gm_distance(arma::mat params1, arma::vec weights1,
                     arma::mat params2, arma::vec weights2);

  void sampleGroupAssignments();

  void computeDistances();

  std::tuple<arma::mat, arma::vec> getParamsAndWeights(int group) {
    arma::mat params(baseClusters[group].size() + hdpClusters.size());
    arma::vec weights(params.n_rows);
    int numBase = baseClusters[group].size();
    for (int h = 0; h < numBase; h++) {
      params.row(h) = baseClusters[group][h].getParam().getParam().t();
      weights(h) = baseClusters[group].size();
    }
    for (int h=0; h < hdpClusters.size(); h++) {
      params.row(numBase + h) = hdpClusters[h].getParam().getParam().t();
      weights(numBase + h) = hdpClusters[h].size(group);
    }
    weights /= arma::sum(weights);
    return std::make_tuple(params, weights);
  }

private:
  using MarginalSemiHdp<Hierarchy>::numGroups;
  using MarginalSemiHdp<Hierarchy>::samplesPerGroup;
  using MarginalSemiHdp<Hierarchy>::data;
  using MarginalSemiHdp<Hierarchy>::c;
  using MarginalSemiHdp<Hierarchy>::datum2cluster;
  using MarginalSemiHdp<Hierarchy>::datum2hierarchy;
  using MarginalSemiHdp<Hierarchy>::baseClusters;
  using MarginalSemiHdp<Hierarchy>::hdpClusters;
  using MarginalSemiHdp<Hierarchy>::engine;
  using MarginalSemiHdp<Hierarchy>::sampleAssignmentsForGroup;
};

#include "metropolised_marginal_semi_hdp_imp.hpp"

#endif //CPPMODEL_METROPOLISED_MARGINAL_SEMI_HDP_HPP
