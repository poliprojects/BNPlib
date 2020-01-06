//
// Created by mario on 25/05/19.
//
#ifndef CPPMODEL_SEMI_HDP_BLOCKED_GIBBS_IMP_HPP
#define CPPMODEL_SEMI_HDP_BLOCKED_GIBBS_IMP_HPP

#include "semi_hdp_blocked_gibbs.hpp"


SemiHdpBlockedGibbs::SemiHdpBlockedGibbs(
    int numGroups, int H): BlockedGibbs(numGroups, H) {}

SemiHdpBlockedGibbs::SemiHdpBlockedGibbs(
    std::vector<arma::vec> &_data, int H): BlockedGibbs(_data, H)  {

  gTildeAtoms = std::make_shared<arma::mat>(numComponents, 2);
  gTildeProbas = std::make_shared<arma::vec>(numComponents);
  numData = 0;
  for (arma::vec dataGroup: data)
    numData += dataGroup.n_elem;

  w_a = 2.0;
  w_b = 2.0;
}

void SemiHdpBlockedGibbs::init() {
  gTildeAtoms->col(0) = stats::rnorm<arma::vec, double>(
      numComponents, 1, 0.0, 10.0);
  gTildeAtoms->col(1) = stats::runif<arma::vec, double>(
      numComponents, 1, 0.5, 1.5);
  *gTildeProbas = stats::rstickbreak(1, gamma, numComponents);
  semiHdpWeight = std::make_shared<double>(0.5);

  clusters.resize(numGroups);

  for (int i = 0; i < numGroups; i++) {
    c(i) = i;
    for (int h = 0; h < numComponents; h++) {
      SemiHdpHierarchy hierarchy(semiHdpWeight, gTildeAtoms,
                                 gTildeProbas);
      hierarchy.samplePrior();
      clusters[i].push_back(Cluster<SemiHdpHierarchy>(hierarchy));
      clusters[i][h].add(data[i][h]);
      datum2cluster[i][h] = h;
    }

    for (int j = numComponents; j < samplesPerGroup[i]; j++) {
      int num = stats::rdiscreteunif(0, numComponents, engine);
      clusters[i][num].add(data[i](j));
      datum2cluster[i][j] = num;
    }

    clusterProbas.row(i) = stats::rstickbreak(
        1.0, alpha, numComponents, engine).t();
  }

}


void SemiHdpBlockedGibbs::maybeMergeClusters() {
  bool merged = false;

  for (int i = 0; i < numGroups; i++) {
    std::unordered_map<int, std::vector<int>> comp2clusteridx;
    for (int h = 0; h < numComponents; h++) {
      if (clusters[i][h].isEmpty())
        continue;

      int component = clusters[i][h].getParam().getGTildeComponent();
      if (component > - 1) {
        comp2clusteridx[component].push_back(h);
      }
    }

    for (auto it: comp2clusteridx) {
      if (it.second.size() < 2)
        continue;

      merged = true;
      int gTildeComp = it.first;
      int h = it.second[0];
      std::unordered_set<int> clustersToMerge;
      Cluster<SemiHdpHierarchy>* firstCluster = &clusters[i][h];
      for (int j=1; j < it.second.size(); j++) {
        int currH = it.second[j];
        clustersToMerge.insert(currH);
        auto data = clusters[i][currH].getData();
        if (! clusters[i][currH].isEmpty() )
          firstCluster->merge(clusters[i][currH]);
        clusters[i][currH].empty();
      }

      arma::uvec groups = arma::find(c == i);
      for (int g : groups) {
        for (int j = 0; j < samplesPerGroup[g]; j++) {
          if (clustersToMerge.find(datum2cluster[g][j]) != clustersToMerge.end()) {
            datum2cluster[g][j] = h;
          }
        }
      }
    }

  }

}

void SemiHdpBlockedGibbs::sampleGtilde() {
  // Atoms
  std::vector<arma::vec> idx2data(numComponents);
  arma::vec idx2count(numComponents, arma::fill::zeros);
  for (int h = 0; h < numComponents; h++) {
    arma::vec curr;
    curr.empty();
    idx2data[h] = curr;
  }

  for (int i = 0; i < numGroups; i++) {
    for (int h = 0; h < numComponents; h++) {
      int component = clusters[i][h].getParam().getGTildeComponent();
      if ((component > -1) && !(clusters[i][h].isEmpty())) {
        idx2data[component] = arma::join_cols(
            idx2data[component], arma::vec(clusters[i][h].getData()));
        idx2count(component) += clusters[i][h].size();
      }
    }
  }

  for (int h = 0; h < numComponents; h++) {
    arma::vec datavec = idx2data[h];
    // FIXME here are hardcoded G_00 hyperparams
    // assert (datavec.n_elem <= data[0].size() + data[1].size());
    arma::vec params = normalGammaUpdate(datavec, 0.0, 2.0, 2, 1.0);
    arma::vec theta(2);
    theta(1) = 1 / stats::rgamma(params(1), 1 / params(2));
    theta(0) = stats::rnorm(params(0), sqrt(1 / theta(1) * params(3)));
    gTildeAtoms->row(h) = theta.t();
  }

  // Probabilities
  *gTildeProbas = _updateClusterProbas(idx2count, gamma);
}

void SemiHdpBlockedGibbs::sampleSemiHdpWeight() {
  int cnt = 0;
  for (int i = 0; i < numGroups; i++)  {
    for (int h = 0; h < numComponents; h++) {
      int component = clusters[i][h].getParam().getGTildeComponent();
      if ((component > -1) && !(clusters[i][h].isEmpty())) {
        cnt += clusters[i][h].size();
      }
    }
  }
  assert(cnt <= numData);
  *semiHdpWeight = stats::rbeta(
      w_a + arma::sum(samplesPerGroup) - cnt, w_b + cnt);
}


void SemiHdpBlockedGibbs::collect(MemoryCollectorByParam *collector) {
  collector->collect("w", *semiHdpWeight);
  collector->collect("gTildeAtoms", *gTildeAtoms);
  collector->collect("gTildeProbas", *gTildeProbas);
  BlockedGibbs::collect(collector);
}

#endif  // CPPMODEL_SEMI_HDP_BLOCKED_GIBBS_IMP_HPP
