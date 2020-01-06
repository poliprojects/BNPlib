//
// Created by mario on 25/05/19.
//

#ifndef CPPMODEL_SEMI_HDP_BLOCKED_GIBBS_HPP
#define CPPMODEL_SEMI_HDP_BLOCKED_GIBBS_HPP

#include "blocked_gibbs.hpp"
#include "semi_hdp_hierarchy.hpp"



class SemiHdpBlockedGibbs: public BlockedGibbs<SemiHdpHierarchy>{
 protected:
  std::shared_ptr<double> semiHdpWeight;
  std::shared_ptr<arma::mat> gTildeAtoms;
  std::shared_ptr<arma::vec> gTildeProbas;
  int numData;

  double w_a; // parameters of w
  double w_b;


 public:

  SemiHdpBlockedGibbs(int numGroups, int H);

  SemiHdpBlockedGibbs(std::vector<arma::vec> & _data, int H);

  void init();

  void sample(bool warmup=false) {
    sampleClusters();
    maybeMergeClusters();
    sampleClusterProbs();
    sampleAssignments();
    if (!warmup)
        sampleGroupAssignments();
    sampleGtilde();
    sampleSemiHdpWeight();
  }

  void maybeMergeClusters();

  void sampleGtilde();

  void sampleSemiHdpWeight();

  void collect(MemoryCollectorByParam* collector);
};

#include "semi_hdp_blocked_gibbs_imp.hpp"

#endif //CPPMODEL_SEMI_HDP_BLOCKED_GIBBS_HPP
