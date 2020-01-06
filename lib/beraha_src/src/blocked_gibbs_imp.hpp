//
// Created by mario on 30/04/19.
//

#ifndef CPPMODEL_BLOCKED_GIBBS_IMP_HPP
#define CPPMODEL_BLOCKED_GIBBS_IMP_HPP

#include "blocked_gibbs.hpp"
#include "options.hpp"

template <typename  Hierarchy>
BlockedGibbs<Hierarchy>::BlockedGibbs(int numGroups, int H):
    numGroups(numGroups), numComponents(H) {

  c.resize(numGroups);
  clusterProbas.resize(numGroups, numComponents);
  clusters.resize(numGroups);

  for (int i = 0; i < numGroups; i++) {
    clusters[i].resize(numComponents);
  }

  // FIXME Fixing HyperParams
  alpha = 1.0;
  gamma = 1.0;
  alpha_0 = 1.0;
}

template <typename  Hierarchy>
BlockedGibbs<Hierarchy>::BlockedGibbs(std::vector<arma::vec> & _data, int H):
    data(_data), numComponents(H) {

  numGroups = data.size();

  c.resize(numGroups);
  samplesPerGroup.set_size(numGroups);
  datum2cluster.resize(numGroups);
  clusterProbas.resize(numGroups, numComponents);
  for (int i=0; i < numGroups; i++) {
    samplesPerGroup(i) = data[i].n_elem;
    datum2cluster[i].resize(data[i].n_elem);
  }

  // FIXME Fixing HyperParams
  alpha = 1.0;
  gamma = 1.0;
  alpha_0 = 1.0;
}

template <typename  Hierarchy>
void BlockedGibbs<Hierarchy>::init() {
  clusters.resize(numGroups);

  for (int i = 0; i < numGroups; i++) {
    c(i) = i;
    for (int h = 0; h < numComponents; h++) {
      NormalConjugateHierarchy hierarchy;
      clusters[i].push_back(Cluster<Hierarchy>(hierarchy));
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

template <typename  Hierarchy>
void BlockedGibbs<Hierarchy>::samplePseudoPrior(int burnin, int numSteps) {
  pseudoPriorCollector.setMaxNumSteps(numSteps + 10);
  pseudoPriorCollector.setNumSteps(numSteps);
  numPseudoPriorIters = numSteps;
  usePseudoPrior = true;

  for (int i = 0; i < burnin; i++)
    sample(true);

  for (int i = 0; i < numSteps; i++) {
    sample(true);
    pseudoPriorCollector.collect("clusterProbas", clusterProbas);

    for (int g = 1; g < numGroups; g++) {
      arma::mat clusterParams(numComponents, 2);
      for (int h = 0; h < numComponents; h++)
        clusterParams.row(h) = clusters[g][h].getParam().getParam().t();

      pseudoPriorCollector.collect(
        "clusterParams" + std::to_string(g), clusterParams);
    }
  }
}


template <typename  Hierarchy>
arma::uvec BlockedGibbs<Hierarchy>::samplePolyaUrn(double concentration, int size) {
  arma::uvec out(size);
  out.zeros();
  for (int i = 1; i < size; i++) {
    arma::vec probs(i+1);
    probs(i) = concentration / (concentration + i);
    for (int j = 0; j < i; j++)
      probs(j) = 1 / (concentration + i);

    int pos = stats::rdiscrete(probs, engine);
    if (pos < i)
      out(i) = out(pos);
    else
      out(i) = out.max() + 1;
  }
  return out;
}

// Useful stuff:
//
// std::vector< std::vector<Cluster<Hierarchy>> > clusters;  (i.e. clusters[i][j])
//
// void sample(bool warmpup=false) {   // for vectors phi, c, p
//   sampleClusters();      // phi
//   sampleClusterProbs();  // p
//   sampleAssignments();   // c
//   if (! warmpup)
//     sampleGroupAssignments();
// }
//
// std::vector<std::vector<int>> datum2cluster;  is the labels vector (the true one, c)
// arma::uvec c;  contains the labels of the GROUPS

template <typename  Hierarchy>
void BlockedGibbs<Hierarchy>::sampleClusters() {
  for (int i=0; i < numGroups; i++) { // for each cluster
    arma::uvec _isUsed = arma::find(c==i);
    bool isUsed = _isUsed.n_elem > 0;
    if (! isUsed && usePseudoPrior) { // _isUsed is empty and ...
      int iternum = stats::rdiscreteunif(0, numPseudoPriorIters-1, engine);
      	// discrete uniform distro with support from 0 to num...-1
      arma::mat probas = pseudoPriorCollector.get("clusterProbas", iternum);
      arma::mat params = pseudoPriorCollector.get(
        "clusterParams" + std::to_string(i), iternum);
      clusterProbas.row(i) = arma::conv_to<arma::vec>::from(probas.row(i)).t();
      for (int h = 0; h < numComponents; h++)
        clusters[i][h].setParam(Hierarchy(params.row(h).t()));
    } else {
      // #pragma omp parallel for
      for (int h = 0; h < numComponents; h++) {
        clusters[i][h].sample();
      }
    }
  }
}

template <typename  Hierarchy>
void BlockedGibbs<Hierarchy>::sampleAssignmentsForGroup(
    int group, bool flush) {

  int component = c(group); // "component"??? It's the label for group g
  for (int j = 0; j < samplesPerGroup(group); j++) { // for each unit in the group
    arma::vec probas(numComponents);
    double datum = data[group](j);

    int oldAssignment = datum2cluster[group][j];
    arma::mat p(numComponents, 2);
    // #pragma omp parallel for
    for (int h = 0; h < numComponents; h++) {
      probas(h) = (clusterProbas(component, h) + 1e-5) * std::exp(
          clusters[component][h].loglike(datum));
      p.row(h) = clusters[component][h].getParam().getParam().t();
    }
    probas /= arma::sum(probas);

    int newAssignment = stats::rdiscrete(probas);
    if (flush) { // if called from sampleGroupAssignments()
      clusters[component][newAssignment].add(datum);
      datum2cluster[group][j] = newAssignment;
    } else if (newAssignment != oldAssignment) { // if called from sampleAssignments()
      clusters[component][oldAssignment].remove(datum);
      clusters[component][newAssignment].add(datum);
      datum2cluster[group][j] = newAssignment;
    }
  }
}


template <typename  Hierarchy>
void BlockedGibbs<Hierarchy>::sampleAssignments() { // uses flush=false
  // #pragma omp parallel for
  for (int i = 0; i < numGroups; i++) {
    sampleAssignmentsForGroup(i, false);
  }
}


template <typename  Hierarchy>
void BlockedGibbs<Hierarchy>::sampleClusterProbs() {
  // #pragma omp parallel for
  for (int i = 0; i < numGroups; i++) {
    arma::vec clusterSizes(numComponents);
    for (int h = 0; h < numComponents; h++)
      clusterSizes(h) = clusters[i][h].size();

    clusterProbas.row(i) = _updateClusterProbas(clusterSizes, alpha).t();
  }

}

template <typename  Hierarchy>
arma::vec BlockedGibbs<Hierarchy>::_updateClusterProbas(arma::vec sizes, double alpha) {
  arma::vec avec = sizes.head(numComponents - 1) + 1;
  avec += 1;

  // b_k = alpha + \sum_{l=k+1}^H M_l
  arma::vec bvec = arma::reverse(arma::cumsum(
      arma::reverse(sizes.tail(numComponents - 1))));
  bvec += alpha;

  return stats::rstickbreak(avec, bvec, engine);
}


template <typename  Hierarchy>
void BlockedGibbs<Hierarchy>::sampleGroupAssignments() { // called if (! warmpup),
  if (numGroups == 1) {                                  // uses flush=true
    return;
  }

  for (int g=1; g < numGroups; g++) {
    int prev = c(g);
    arma::vec logprobas(numGroups, arma::fill::zeros);
    // TODO ASSIGN PRIOR!
    // logprobas(0) = - log(1 + alpha_0);
    // logprobas(1) = log(alpha_0) - log(1 + alpha_0);
    int numSamples = samplesPerGroup(1);
    for (int i = 0; i < numGroups; i++) {
      for (int j = 0; j < numSamples; j++) {
        logprobas(i) += logLikeForComponent(data[1](j), i);
      }
    }
    arma::vec probas = arma::exp(logprobas);
    probas /= arma::sum(probas);
    c(g) = stats::rdiscrete(probas, engine);

    if (c(g) != prev) { // if a new LABEL for group g was extracted
      for (int j = 0; j < samplesPerGroup(g); j++) {
        double datum = data[g](j);
        clusters[prev][datum2cluster[g][j]].remove(datum);
        datum2cluster[g][j] = -1;
        	// you're destroying the whole group g
      }
      sampleAssignmentsForGroup(g, true);
    }
  }
}


template <typename  Hierarchy>
double BlockedGibbs<Hierarchy>::logLikeForComponent(double x, int component) {
  arma::vec logprobas(numComponents);
  // #pragma omp parallel for
  for (int h=0; h < numComponents; h++) {
    logprobas(h) = log(clusterProbas(component, h)) + \
                   clusters[component][h].loglike(x);
  }
  return logSumExp(logprobas);

}


template <typename  Hierarchy>
void BlockedGibbs<Hierarchy>::collect(MemoryCollectorByParam *collector) {
  collector->collect("c", arma::conv_to<arma::vec>::from(c));
  for (int i = 0; i < numGroups; i++) {
    collector->collect(
        "s."+std::to_string(i+1), stdToArma(datum2cluster[i]));
  }
  collector->collect("clusterProbas", clusterProbas);
  for (int i = 0; i < numGroups; i++) {
    std::string name = "cluster." + std::to_string(i + 1) + ".";
    for (int h = 0; h < numComponents; h++) {
      collector->collect(
          name + std::to_string(h),
          clusters[i][h].getParam().getParam());
    }
  }
}


template <typename Hierarchy>
arma::mat BlockedGibbs<Hierarchy>::predictiveSamples(
    const MemoryCollectorByParam& collector) {
  int numSteps = collector.getNumSteps();
  arma::mat out(numSteps, numGroups);

  for (int iternum = 0; iternum < numSteps; iternum++) {
    c = arma::conv_to<arma::uvec>::from(collector.get(std::string("c"), iternum));

    // cluster probas
    clusterProbas = collector.get(std::string("clusterProbas"), iternum);

    // clusters
    for (int i = 0; i < numGroups; i++) {
      for (int h = 0; h < numComponents; h++) {
        arma::vec param = collector.get(
            "cluster." + std::to_string(i + 1) + "." + std::to_string(h),
            iternum);
        clusters[i][h].setParam(Hierarchy(param));
      }
    }

    for (int i = 0; i < numGroups; i++) {
      int comp = stats::rdiscrete(clusterProbas.row(c(i)), engine);
      out(iternum, i) = clusters[c(i)][comp].getParam().predict();
    }
  }
  return out;
}


template <typename  Hierarchy>
void BlockedGibbs<Hierarchy>::print() {
  for (int i=0; i < numGroups; i++) {
    std::cout << "*** Group: " << i << " ***" << std::endl;
    for (int h = 0; h < numComponents; h++) {
      arma::vec data(clusters[i][h].getData());
      if (data.n_elem > 0) {
        std::cout << "Cluster # " << h << std::endl;
        std::cout << "Parameter: "; clusters[i][h].getParam().getParam().t().print();
        std::cout << "Proba: " << clusterProbas(i, h) << std::endl;
        // std::cout << "Data: "; data.t().print();
      }
    }
  }
}


template <typename  Hierarchy>
int BlockedGibbs<Hierarchy>::getNumGroups() const {
  return numGroups;
}


template <typename  Hierarchy>
const arma::uvec &BlockedGibbs<Hierarchy>::getSamplesPerGroup() const {
  return samplesPerGroup;
}


template <typename  Hierarchy>
const arma::uvec &BlockedGibbs<Hierarchy>::getC() const {
  return c;
}


template <typename  Hierarchy>
const arma::mat &BlockedGibbs<Hierarchy>::getClusterProbas() const {
  return clusterProbas;
}

#endif // CPPMODEL_BLOCKED_GIBBS_IMP_HPP
