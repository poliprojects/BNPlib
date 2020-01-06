//
// Created by mario on 27/06/19.
//

#ifndef CPPMODEL_MARGINAL_SEMI_HDP_IMP_HPP
#define CPPMODEL_MARGINAL_SEMI_HDP_IMP_HPP

#include "marginal_semi_hdp.hpp"

template <typename Hierarchy>
MarginalSemiHdp<Hierarchy>::MarginalSemiHdp(
    std::vector<arma::vec> & _data): data(_data) {
  numGroups = data.size();
  std::cout << "Numgroups: " << numGroups << std::endl;
  samplesPerGroup.resize(numGroups);
  datum2cluster.resize(numGroups);
  datum2hierarchy.resize(numGroups);
  c.resize(numGroups);
  for (int i = 0; i < numGroups; i++) {
    samplesPerGroup(i) = data[i].size();
    c(i) = i;
    datum2cluster[i].resize(samplesPerGroup(i));
    datum2hierarchy[i].resize(samplesPerGroup(i));
  }
  nDataToHdpFromRestaurant = arma::uvec(numGroups, arma::fill::zeros);
  alpha = 2.0;
  gamma = 2.0;
  alpha_0 = 1.0;
  w_a = 2.0;
  w_b = 2.0;
  etas = arma::vec(numGroups, arma::fill::ones);
  omegas = stats::rdirichlet(etas, engine);
}


template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::init() {
  int CLUSTER_RATIO = 10;
  int HDP_CLUSTER_RATIO = 10;
  baseClusters.resize(numGroups);
  pseudoBaseClusters.resize(numGroups);
  int numHdpClus = arma::sum(samplesPerGroup) / (HDP_CLUSTER_RATIO + numGroups);
  semiHdpWeight = 0.5;
  for (int h0 = 0; h0 < numHdpClus; h0++) {
    Hierarchy hierarchy(G00Params);
    hdpClusters.push_back(Cluster<Hierarchy>(hierarchy, numGroups));
  }

  // Initialize clusters
  for (int i = 0; i < numGroups; i++) {
    int numClus = samplesPerGroup(i) / CLUSTER_RATIO + \
                  (samplesPerGroup(i) % CLUSTER_RATIO != 0);
    int h;
    for (h = 0; h < numClus; h++) {
      Hierarchy hierarchy(G0Params);
      baseClusters[i].push_back(Cluster<Hierarchy>(hierarchy, numGroups));
      baseClusters[i][h].add(data[i][h], i);
      datum2cluster[i][h] = h;
      datum2hierarchy[i][h] = false;
    }

    for (h = numClus; h < samplesPerGroup(i); h++) {
      if (stats::rbern(0.5, engine)) {
        int comp = stats::rdiscreteunif(0, numHdpClus, engine);
        hdpClusters[comp].add(data[i][h], i);
        datum2cluster[i][h] = comp;
        datum2hierarchy[i][h] = true;
        nDataToHdpFromRestaurant(i) += 1;
      } else {
        int comp = stats::rdiscreteunif(0, numClus, engine);
        baseClusters[i][comp].add(data[i][h], i);
        datum2cluster[i][h] = comp;
        datum2hierarchy[i][h] = false;
      }
    }
  }
  // TODO initialize numTables, betas
  nDataToHdp = arma::sum(nDataToHdpFromRestaurant);
  numTables = arma::vec(hdpClusters.size(), arma::fill::ones);
  relabel();
  arma::vec alpha_param = arma::vec(hdpClusters.size(), arma::fill::ones);
  alpha_param = arma::join_cols(alpha_param, arma::vec{gamma});
  betas = arma::vec(hdpClusters.size() + 1, arma::fill::ones);
  betas /= arma::sum(betas);
  std::cout << "Betas init: "; betas.t().print();
}


template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::samplePseudoPrior(
    int burnin, int numSteps) {
  // TODO take care of numTables and betas
  pseudoPriorCollector.setMaxNumSteps(numSteps + 10);
  pseudoPriorCollector.setNumSteps(numSteps);
  numPseudoPriorIters = numSteps;
  usePseudoPrior = true;

  for (int i = 0; i < burnin; i++)
    sample(true);

  for (int i = 0; i < numSteps; i++) {
    sample(true);
    pseudoPriorCollector.newIter();

    for (int r=1; r < numGroups; r++) {
      arma::mat baseClusterParams(baseClusters[r].size(), 2);
      arma::vec numCustomersPerCluster(baseClusters[r].size());
      for (int h = 0; h < baseClusters[r].size(); h++) {
        baseClusterParams.row(h) = baseClusters[r][h].getParam().getParam().t();
        numCustomersPerCluster(h) = baseClusters[r][h].size();
      }
      pseudoPriorCollector.collect(
        "baseClusterParams." + std::to_string(r), baseClusterParams);
      pseudoPriorCollector.collect(
          "numCustomersPerCluster." + std::to_string(r), numCustomersPerCluster);
    }

    arma::mat hdpClusterParams(hdpClusters.size(), 2);
    arma::mat numCustomerPerHdpCluster(hdpClusters.size(), numGroups);
    for (int h = 0; h < hdpClusters.size(); h++) {
      hdpClusterParams.row(h) = hdpClusters[h].getParam().getParam().t();
      numCustomerPerHdpCluster.row(h) = arma::conv_to<arma::rowvec>::from(
        hdpClusters[h].getNumDataFromRestaurant());
    }
    pseudoPriorCollector.collect("hdpClusterParams", hdpClusterParams);
    pseudoPriorCollector.collect("numCustomerPerHdpCluster", numCustomerPerHdpCluster);
  }
}


template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::sampleClusters() {
  for (auto cc: hdpClusters)
    cc.sample();

  #pragma omp parallel for
  for (int i=0; i < numGroups; i++) {
    arma::uvec _isUsed = arma::find(c==i);
    bool isUsed = _isUsed.n_elem > 0;
    if (! isUsed && usePseudoPrior) {
      // TODO take care of numTables and betas
      pseudoBaseClusters[i].resize(0);
      pseudoHdpClusters.resize(0);

      int iternum = stats::rdiscreteunif(0, numPseudoPriorIters-1, engine);
      arma::mat baseClusterParams = pseudoPriorCollector.get(
          "baseClusterParams." + std::to_string(i), iternum);
      arma::mat numCustomersPerCluster = pseudoPriorCollector.get(
          "numCustomersPerCluster." + std::to_string(i), iternum);
      pseudoBaseClusters.reserve(baseClusterParams.n_rows);
      for (int h=0; h < baseClusterParams.n_rows; h++) {
        Hierarchy param(baseClusterParams.row(h).t(),
                        G00Params);
        Cluster<Hierarchy> clus(param, numGroups);
        clus.setNumPseudoElems(numCustomersPerCluster(h));
        pseudoBaseClusters[i].push_back(clus);
      }
      baseClusters[i] = pseudoBaseClusters[i];

      // TODO change here if we want to add also HdpClusters
      // coming from the pseudoprior
//      arma::mat hdpClusterParams = pseudoPriorCollector.get(
//          "hdpClusterParams", iternum);
//      arma::mat numCustomerPerHdpCluster = pseudoPriorCollector.get(
//          "numCustomerPerHdpCluster", iternum);
//      pseudoHdpClusters.reserve(hdpClusterParams.n_rows);
//      for (int h=0; h < hdpClusterParams.n_rows; h++) {
//        Hierarchy param(hdpClusterParams.row(h).t());
//        Cluster<Hierarchy> clus(param, numGroups);
//        arma::uvec ndata = arma::conv_to<arma::uvec>::from(
//            numCustomerPerHdpCluster.row(h).t());
//        clus.setNPseudoDataFromRestaurant(ndata);
//        pseudoHdpClusters.push_back(clus);
//      }
//
//      for (auto clus: pseudoHdpClusters)
//        hdpClusters.push_back(clus);
//
//      numTables = arma::join_cols(
//        numTables, arma::vec(pseudoHdpClusters.size(), arma::fill::ones));
//      sampleHdpLatentBetas();
    } else {
    for (auto cc: baseClusters[i])
      cc.sample();
    }
  }
  // TODO check here that we are not sampling from the
  // prior and removing all the good stuff
}


template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::sampleAssignmentsForGroup(
    int group, bool flush) {
  int component = c(group);

  int nDataToBase = 0;
  for (auto cc: baseClusters[component])
    nDataToBase += cc.size();

  int newAssignment, oldAssignment;

  for (int j=0; j < samplesPerGroup(group); j++) {
    oldAssignment = datum2cluster[group][j];
    double datum = data[group][j];
    if (datum2hierarchy[group][j]) {
      newAssignment = _assignDatumToHdp(datum, group, component, true);
      datum2cluster[group][j] = newAssignment;
      if (oldAssignment != -1) {
        hdpClusters[oldAssignment].remove(datum, component);
      }
    } else {
      // CRP Withing group
      newAssignment = _assignDatumToRestaurant(datum, group, component, true);
      datum2cluster[group][j] = newAssignment;
      if (oldAssignment != -1)
        baseClusters[component][oldAssignment].remove(datum, component);
    }
    assert(betas.n_elem == hdpClusters.size() + 1);
  }
}


template <typename Hierarchy>
int MarginalSemiHdp<Hierarchy>::_assignDatumToRestaurant(
    double datum, int group, int restaurant, bool addNew) {
  int nclus = baseClusters[restaurant].size();
  arma::vec probas(nclus + 1, arma::fill::zeros);
  for (int h = 0; h < nclus; h++) {
    probas[h] = exp(baseClusters[restaurant][h].logProbaForDatum(datum));
  }
  if (addNew) {
    probas[nclus] = alpha * exp(Hierarchy::marginalLogLike(datum, G0Params));
  }
  probas /= arma::sum(probas);
  int newAssignment = stats::rdiscrete(probas, engine);

  if (newAssignment == baseClusters[restaurant].size()) {
    baseClusters[restaurant].push_back(Cluster<Hierarchy>(numGroups));
    baseClusters[restaurant][newAssignment].add(datum, restaurant, true);
  } else {
    baseClusters[restaurant][newAssignment].add(datum, restaurant);
  }
  return newAssignment;
}


template <typename Hierarchy>
int MarginalSemiHdp<Hierarchy>::_assignDatumToHdp(
    double datum, int group, int restaurant, bool addNew) {
  int nclus = hdpClusters.size();
  arma::vec probas(nclus + 1, arma::fill::zeros);
  for (int h = 0; h < nclus; h++) {
    probas(h) = _probaForHdpClus(h, datum, restaurant);
  }
  if (addNew)
    probas(nclus) = _probaForNewHdpClus(datum, restaurant);

  probas /= arma::sum(probas);
  int newAssignment = stats::rdiscrete(probas, engine);

  if (newAssignment == hdpClusters.size()) {
    hdpClusters.push_back(Cluster<Hierarchy>(numGroups));
    hdpClusters[newAssignment].add(datum, restaurant, true);
    numTables = arma::join_cols(numTables, arma::vec{1});
    sampleHdpLatentBetas();
  } else {
    hdpClusters[newAssignment].add(datum, restaurant);
  }
  assert(betas.n_elem == hdpClusters.size() + 1);
  return newAssignment;
}


template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::sampleAssignments() {
  // First we sample numTables and beta
  nDataToHdpFromRestaurant = arma::uvec(numGroups, arma::fill::zeros);
  for (auto cc: hdpClusters)
    nDataToHdpFromRestaurant += cc.getNumDataFromRestaurant();

  sampleHdpLatentNtables();
  sampleHdpLatentBetas();

  for (int i = 0; i < numGroups; i++) {
    sampleAssignmentsForGroup(i, false);
  }
}

template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::sampleHdpLatentNtables() {
  numTables = arma::vec(hdpClusters.size(), arma::fill::ones);
  for (int i=0; i < numGroups; i++) {
    arma::vec curr(hdpClusters.size(), arma::fill::zeros);
    for (int h=0; h < hdpClusters.size(); h++) {
      int numCustomers = hdpClusters[h].size(i);
      double beta_h = betas(h);
      arma::vec probas(numCustomers, arma::fill::zeros);
      for (int m=0; m < probas.n_elem; m++) {
        double s = 1.0 * stirling(numCustomers, m + 1);
        double gammas =  exp(
            lgamma(alpha * beta_h) - lgamma(alpha * beta_h + numCustomers));
        double power = std::pow(alpha * betas[h], m + 1);
        double proba = s * gammas * power;
        probas(m) = proba;
      }
      if (arma::sum(probas) > 0) {
        probas /= arma::sum(probas);
        curr(h) = stats::rdiscrete(probas, engine);
      } else {
        // This happens because we might have empty clusters!
        curr(h) = 0;
      }
    }
    numTables += curr;
  }
}

template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::sampleHdpLatentBetas() {
  betas = stats::rdirichlet(arma::join_cols(numTables, arma::vec{gamma}));
}

template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::sampleHieararchy() {
  for (arma::uword i=0; i < numGroups; i++) {
    int nDataToBase = 0;
    for (auto cc: baseClusters[c(i)])
      nDataToBase += cc.size();

    int nDataToHdp = 0;
    for (auto cc: hdpClusters)
      nDataToHdp += cc.size((int) c(i));

    std::vector<bool> currs(samplesPerGroup(i));
    #pragma omp parallel for
    for (int j=0; j < samplesPerGroup(i); j++) {
      double datum = data[i][j];
      arma::vec probas(2, arma::fill::zeros);

      arma::vec prob0(baseClusters[c(i)].size() + 1);
      for (int h=0; h < baseClusters[c(i)].size(); h++)
        prob0(h) = baseClusters[c(i)][h].logProbaForDatum(datum);

      prob0(baseClusters[c(i)].size()) = log(alpha) + \
        Hierarchy::marginalLogLike(datum, G0Params);
      prob0 -= log(nDataToBase + alpha);
      probas(0) = logSumExp(prob0);

      arma::vec prob1(hdpClusters.size() + 1);
      for (int h=0; h < hdpClusters.size(); h++)
        prob1(h) = log(_probaForHdpClus(h, datum, c(i)));

      prob1(hdpClusters.size()) = log(_probaForNewHdpClus(datum, c(i)));
      prob1 -= log(nDataToHdp + alpha);
      probas(1) = logSumExp(prob1);
      probas = arma::exp(probas);
      probas /= arma::sum(probas);

      // probas[1] = 0.0;
      currs[j] = (stats::rbern(probas[1], engine) == 1);
    }
    for (int j=0; j < samplesPerGroup(i); j++) {
      bool curr = currs[j];
      bool prev = datum2hierarchy[i][j];
      double datum = data[i][j];
      if (prev && (!curr)) {
        // We are removing one datum from the hdp tables, where do we
        // assign it ?
        int oldAssignment = datum2cluster[i][j];
        hdpClusters[oldAssignment].remove(datum, c(i));
        int newAssignment = _assignDatumToRestaurant(datum, i, c(i), true);
        datum2cluster[i][j] = newAssignment;
      } else if (curr && !(prev)) {
        // We are adding one datum to "some" hdp table, which one?
        int oldAssignment = datum2cluster[i][j];
        baseClusters[c(i)][oldAssignment].remove(datum, c(i));
        int newAssignment = _assignDatumToHdp(datum, i, c(i), false);
        datum2cluster[i][j] = newAssignment;
      }
      datum2hierarchy[i][j] = curr;
    }
  }
}

template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::relabel() {
  // TODO take care of betas, numTables
  // HDP
  std::vector<int> toRemoveHdp;
  for (int h = 0; h < hdpClusters.size(); h++) {
    if (hdpClusters[h].isEmpty())
      toRemoveHdp.push_back(h);
  }
  for (auto index=toRemoveHdp.rbegin(); index != toRemoveHdp.rend(); index ++) {
    int ind = *index;
    hdpClusters.erase(hdpClusters.begin() + ind);
    numTables.shed_row(ind);
  }

  sampleHdpLatentBetas();

  // BASE
  for (int g=0; g < numGroups; g++) {
    std::vector<int> toRemove;
    for (int h=0; h < baseClusters[g].size(); h++) {
      if (baseClusters[g][h].isEmpty()) {
        toRemove.push_back(h);
      }
    }
    for (auto index=toRemove.rbegin(); index != toRemove.rend(); index ++) {
      int ind = *index;
      baseClusters[g].erase(baseClusters[g].begin() + ind);
    }
    arma::uvec dataGroups = arma::find(c == g);
    for (int i = 0; i < dataGroups.n_elem; i++) {
      int group = dataGroups[i];
      for (int j=0; j < samplesPerGroup[group]; j++) {
        if (!datum2hierarchy[group][j]) {
          datum2cluster[group][j] -= std::count_if(
              toRemove.begin(), toRemove.end(),
              [this, &group, &j](int k) {return k < datum2cluster[group][j];});
        } else {
          datum2cluster[group][j] -= std::count_if(
              toRemoveHdp.begin(), toRemoveHdp.end(),
              [this, &group, &j](int k) {return k < datum2cluster[group][j];});
        }
      }
    }
  }
}

template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::sampleGroupAssignments() {
  if (numGroups == 1) {
    return;
  }

  for (int g=1; g < numGroups; g++) {
    int prev = c(g);
    arma::vec logprobas = arma::log(omegas);
    int numSamples = samplesPerGroup(g);
    for (int j = 0; j < numSamples; j++) {
      arma::vec curr(numGroups);
      #pragma omp parallel for
      for (int i = 0; i < numGroups; i++) {
        curr(i) = logLikeForComponent(data[g][j], i);
      }
      logprobas += curr;
    }

    arma::vec probas = arma::exp(logprobas);

    probas /= arma::sum(probas);
    c(g) = stats::rdiscrete(probas, engine);

    if (c(g) != prev) {
      // std::cout << "Old: " << prev << ", New: " <<
      //   c(g) << ", logProbas: "; logprobas.t().print();
      for (int j = 0; j < samplesPerGroup(g); j++) {
        double datum = data[g](j);
        if (datum2hierarchy[g][j])
          hdpClusters[datum2cluster[g][j]].remove(datum, prev);
        else
          baseClusters[prev][datum2cluster[g][j]].remove(datum, prev);
        datum2cluster[g][j] = -1;
      }
      sampleAssignmentsForGroup(g, true);
  //    relabel();
    }
  }
}


template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::sampleSemiHdpWeight() {
  nDataToHdp = 0;
  for (int i=0; i < numGroups; i++)
    for (bool z : datum2hierarchy[i])
      nDataToHdp += z;

  semiHdpWeight = stats::rbeta(
      w_a + arma::sum(samplesPerGroup) - nDataToHdp, w_b + nDataToHdp);
}


template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::sampleHypers() {
  int numDistinctClusters = 0;
  for (auto rest: baseClusters)
    numDistinctClusters += rest.size();
  arma::mat uniqueRestClusters(numDistinctClusters, 2);
  int row = 0;
  for (auto rest: baseClusters) {
    for (auto cc: rest) {
      uniqueRestClusters.row(row) = cc.getParam().getParam().t();
      row += 1;
    }
  }
  G0Params.update(uniqueRestClusters);

  arma::mat uniqueHdpClusters(hdpClusters.size(), 2);
  for (int row=0; row < hdpClusters.size(); row++) {
    uniqueHdpClusters.row(row) = hdpClusters[row].getParam().getParam().t();
  }

  G00Params.update(uniqueHdpClusters);

  // now update omegas
  arma::vec cnts(numGroups, arma::fill::zeros);
  for (int i=0; i < numGroups; i++) {
    cnts(i) = arma::sum(c == i);
  }
  omegas = stats::rdirichlet(etas + cnts, engine);
}


template <typename Hierarchy>
double MarginalSemiHdp<Hierarchy>::logLikeForComponent(
    double x, int component) {
  arma::vec probas(baseClusters[component].size() + hdpClusters.size(),
                   arma::fill::zeros);
  int nDataFromRest = 0;
  for (auto cc: baseClusters[component])
    nDataFromRest += cc.size(true);

  for (auto cc: hdpClusters)
    nDataFromRest += cc.size(component);

  int Hr = baseClusters[component].size();
  int H0 = hdpClusters.size();
  // First: component specific clusters
  for (int h = 0; h < Hr; h++) {
    probas(h) = baseClusters[component][h].logProbaForDatum(x);
  }

  // Then: shared clusters
  for (int h=0; h < H0; h++) {
    probas(h + Hr) = log(_probaForHdpClus(h, x, component));
  }

  probas -= log(nDataFromRest + alpha);
  double out = logSumExp(probas);
  return out;
}


//template <typename Hierarchy>
//double MarginalSemiHdp<Hierarchy>::logLikeForComponent(
//    double x, int component) {
//  arma::vec baseLogProbas(baseClusters[component].size() + 1, arma::fill::zeros);
//  arma::vec hdpLogProbas(hdpClusters.size() + 1, arma::fill::zeros);
//  double baseLogProba = 0;
//  double hdpLogProba = 0;
//  int nDataToBase = 0;
//  for (auto cc: baseClusters[component])
//    nDataToBase += cc.size(true);
//
//  // First: component specific clusters
//  for (int h = 0; h < baseClusters[component].size(); h++) {
//    baseLogProbas(h) = baseClusters[component][h].logProbaForDatum(x) - \
//                       log(nDataToBase + alpha - 1);
//  }
//
//  Hierarchy param;
//  baseLogProbas(baseClusters[component].size()) = \
//      log(alpha) + param.marginalLogLike(x) \
//      - log(nDataToBase + alpha - 1);
//
//  baseLogProba = logSumExp(baseLogProbas);
//  // baseLogProba += log(semiHdpWeight + 1e-6);
//
//  // Then: shared clusters
//  for (int h=0; h<hdpClusters.size(); h++) {
//    hdpLogProbas(h) = log(_probaForHdpClus(h, x, component));
//  }
//  hdpLogProbas(hdpClusters.size()) = log(_probaForNewHdpClus(x, component));
//  hdpLogProbas -= arma::sum(nDataToHdpFromRestaurant) -1 + alpha;
//
//  hdpLogProba = logSumExp(hdpLogProbas);
//  // hdpLogProba += log(1 - semiHdpWeight + 1e-6);
//  double out = logSumExp(arma::vec{baseLogProba, hdpLogProba});
//
//  return out;
//}


template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::print() {
  std::cout << "Number of shared clusters: " << hdpClusters.size() << std::endl;
  for (int h = 0; h < hdpClusters.size(); h++) {
    arma::vec data(hdpClusters[h].getData());
    if (data.n_elem > 0) {
      std::cout << "Shared cluster # " << h << std::endl;
      std::cout << "Parameter: ";
      hdpClusters[h].getParam().getParam().t().print();
      std::cout << "Data From Restaurant";
      hdpClusters[h].getNumDataFromRestaurant().t().print();
      std::cout << "Data: ";
      data.t().print();
    }
  }

  for (int i=0; i < numGroups; i++) {
    std::cout << "*** Group: " << i << " ***" << std::endl;
    for (int h = 0; h < baseClusters[i].size(); h++) {
      arma::vec data(baseClusters[i][h].getData());
      if (data.n_elem > 0) {
        std::cout << "Cluster # " << h << std::endl;
        std::cout << "Parameter: ";
        baseClusters[i][h].getParam().getParam().t().print();

        std::cout << "Data: ";
        data.t().print();
      }
    }
  }
}

template <typename Hierarchy>
arma::mat MarginalSemiHdp<Hierarchy>::predictiveSamples(
    const MemoryCollectorByIter& collector) {
  int numSteps = collector.getNumSteps();
  arma::mat out(numSteps, numGroups);
  for (int iternum=0; iternum < numSteps; iternum++) {
    _restoreState(collector, iternum);
    arma::vec curr(numGroups);
    for (int i=0; i < numGroups; i++) {
      bool fromBase = stats::rbern(semiHdpWeight, engine);
      int comp = c(i);
      if (fromBase) {
        arma::vec probas(baseClusters[comp].size() + 1);
        for (int h=0; h < baseClusters[comp].size(); h++) {
          probas(h) = baseClusters[comp][h].size();
        }
        probas(baseClusters[comp].size()) = alpha;
        probas /= arma::sum(probas);
        int clus = stats::rdiscrete(probas, engine);
        if (clus == baseClusters[comp].size()) {
          // std::cout << "Predicting from prior" << std::endl;
          Hierarchy param(G00Params);
          curr(i) = param.predict();
        } else {
          // std::cout << "Predicting from ";
          // baseClusters[comp][clus].getParam().getParam().t().print();
          curr(i) = baseClusters[comp][clus].getParam().predict();
        }
      } else {
        arma::vec probas(hdpClusters.size() + 1);
        for (int h=0; h < hdpClusters.size(); h++) {
          int n0lh = hdpClusters[h].size(comp);
          int n0h = hdpClusters.size();
          probas[h] = 1.0 * n0lh + alpha * n0h / (nDataToHdp + alpha_0);
        }
        probas(hdpClusters.size()) = alpha * alpha_0 / (nDataToHdp + alpha_0);
        probas /= arma::sum(probas);
        int clus = stats::rdiscrete(probas, engine);
        if (clus == hdpClusters.size()) {
          // std::cout << "Predicting from prior" << std::endl;
          Hierarchy param(G00Params);
          curr(i) = param.predict();
        } else {
          // std::cout << "Group: " << i << ", Predicting from ";
          // hdpClusters[clus].getParam().getParam().t().print();
          curr(i) = hdpClusters[clus].getParam().predict();
        }
      }
    }
    out.row(iternum) = curr.t();
  }
  return out;
}

template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::collect(
    MemoryCollectorByIter *collector) {
  collector->collect("c", arma::conv_to<arma::vec>::from(c));

  for (int i = 0; i < numGroups; i++) {
    collector->collect(
        "s."+std::to_string(i+1), stdToArma(datum2cluster[i]));
    collector->collect(
        "z."+std::to_string(i+1), stdToArma(datum2hierarchy[i]));

    arma::mat clusterParamsMat(baseClusters[i].size(), 2);
    arma::vec numClientsMat(baseClusters[i].size(), arma::fill::zeros);
    for (int h = 0; h < baseClusters[i].size(); h++) {
      clusterParamsMat.row(h) = baseClusters[i][h].getParam().getParam().t();
      numClientsMat(h) = baseClusters[i][h].size();
    }
    collector->collect("baseClustersParams." + std::to_string(i + 1), clusterParamsMat);
    collector->collect("numCustomerPerCluster." + std::to_string(i + 1), numClientsMat);
  }

  arma::mat hdpClusterParamsMat(hdpClusters.size(), 2);
  arma::mat numCustomerPerHdpCluster(hdpClusters.size(), numGroups);
  for (int h=0; h < hdpClusters.size(); h++) {
    hdpClusterParamsMat.row(h) = hdpClusters[h].getParam().getParam().t();
    numCustomerPerHdpCluster.row(h) = arma::conv_to<arma::rowvec>::from(
        hdpClusters[h].getNumDataFromRestaurant());
  }
  collector->collect("hdpClusterParams", hdpClusterParamsMat);
  collector->collect("numCustomerPerHdpCluster", numCustomerPerHdpCluster);
  collector->collect("w", arma::vec{semiHdpWeight});

  collector->collect("G0.mu0", arma::vec{G0Params.getMu0()});
  collector->collect("G0.a", arma::vec{G0Params.getA()});
  collector->collect("G0.b", arma::vec{G0Params.getB()});
  collector->collect("G0.lambda", arma::vec{G0Params.getLambda()});

  collector->collect("G00.mu0", arma::vec{G00Params.getMu0()});
  collector->collect("G00.a", arma::vec{G00Params.getA()});
  collector->collect("G00.b", arma::vec{G00Params.getB()});
  collector->collect("G00.lambda", arma::vec{G00Params.getLambda()});
}


template <typename Hierarchy>
void MarginalSemiHdp<Hierarchy>::_restoreState(
    const MemoryCollectorByIter &collector, int iternum) {
  c = arma::conv_to<arma::uvec>::from(collector.get(std::string("c"), iternum));
  semiHdpWeight = arma::conv_to<arma::vec>::from(
    collector.get(std::string("w"), iternum))(0);
  for (int i=0; i < numGroups; i++) {
    baseClusters[i].resize(0);
    arma::mat baseClusterParams = collector.get(
        "baseClustersParams." + std::to_string(i + 1), iternum);
    arma::mat numCustomersPerCluster = collector.get(
        "numCustomerPerCluster." + std::to_string(i + 1), iternum);
    baseClusters[i].reserve(baseClusterParams.n_rows);

    for (int h=0; h < baseClusterParams.n_rows; h++) {
      Hierarchy param(baseClusterParams.row(h).t(),
                      G00Params);
      Cluster<Hierarchy> clus(param, numGroups);
      clus.setNumElems(numCustomersPerCluster(h));
      baseClusters[i].push_back(clus);
    }
  }

  hdpClusters.resize(0);
  arma::mat hdpClusterParams = collector.get(
      "hdpClusterParams", iternum);
  arma::umat numCustomerPerHdpCluster = arma::conv_to<arma::umat>::from(
      collector.get("numCustomerPerHdpCluster", iternum));
  nDataToHdp = arma::accu(numCustomerPerHdpCluster);
  hdpClusters.reserve(hdpClusterParams.n_rows);
  for (int h=0; h < hdpClusterParams.n_rows; h++) {
    Hierarchy param(hdpClusterParams.row(h).t(),
                    G00Params);
    Cluster<Hierarchy> clus(param, numGroups);
    arma::uvec ndata = numCustomerPerHdpCluster.row(h).t();
    clus.setNDataFromRestaurant(ndata);
    hdpClusters.push_back(clus);
  }
}

template <typename Hierarchy>
std::vector<arma::mat> MarginalSemiHdp<Hierarchy>::estimateDensities(
    const std::vector<arma::vec>& grids,
    const MemoryCollectorByIter& collector) {

  std::vector<arma::mat> out(numGroups);
  for (int g=0; g < numGroups; g++) {
    out[g].resize(collector.getNumSteps(), grids[g].n_elem);
  }

  for (int iternum=0; iternum < collector.getNumSteps(); iternum++) {
    _restoreState(collector, iternum);
    for (int g=0; g < numGroups; g++) {
      int comp = c(g);
      int Hr = baseClusters[comp].size();
      int H0 = hdpClusters.size();
      for (int j=0; j < grids[g].size(); j++) {
        double x = grids[g](j);
        out[g](iternum, j) = 0.0;
        arma::vec baseProbas(baseClusters[comp].size() + 1);
        for (int h=0; h < Hr; h++) {
          baseProbas(h) = baseClusters[comp][h].size();
        }
        baseProbas(Hr) = alpha;
        baseProbas /= arma::sum(baseProbas);
        baseProbas = arma::log(baseProbas) + std::log(semiHdpWeight);
        for (int h=0;  h < Hr; h++)
          baseProbas(h) += baseClusters[comp][h].getParam().loglike(x);

        baseProbas(Hr) += Hierarchy::marginalLogLike(x, G0Params);

        arma::vec hdpProbas(hdpClusters.size() + 1);
        for (int h=0; h < H0; h++)
          hdpProbas(h) = hdpClusters[h].size(comp);
        hdpProbas[H0] = gamma;
        hdpProbas /= arma::sum(hdpProbas);
        hdpProbas = arma::log(hdpProbas) + std::log(1 - semiHdpWeight);
        for (int h=0;  h < H0; h++)
          hdpProbas(h) += hdpClusters[h].getParam().loglike(x);

        hdpProbas(H0) += Hierarchy::marginalLogLike(x, G00Params);
        out[g](iternum, j) = std::exp(
          logSumExp(arma::join_cols(baseProbas, hdpProbas)));
      }
    }
  }
  return out;
}

template <typename Hierarchy>
const arma::uvec &MarginalSemiHdp<Hierarchy>::getC() const {
  return c;
}

#endif //CPPMODEL_MARGINAL_SEMI_HDP_IMP_HPP
