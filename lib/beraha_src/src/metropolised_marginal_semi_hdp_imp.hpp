#ifndef CPPMODEL_METROPOLISED_MARGINAL_SEMI_HDP_IMP_HPP
#define CPPMODEL_METROPOLISED_MARGINAL_SEMI_HDP_IMP_HPP

#include "metropolised_marginal_semi_hdp.hpp"

template<typename Hierarchy>
double MetropolisedMarginalSemiHdp<Hierarchy>::gm_distance(
    arma::mat params1, arma::vec weights1,
    arma::mat params2, arma::vec weights2) {
  int N1 = weights1.n_elem;
  int N2 = weights2.n_elem;
  double out;
  for (int i=0; i < N1; i++) {
    double alpha_i = weights1(i);
    double m1 = params1(i, 0);
    double sd1 = params1(i, 1);
    double var1 = sd1 * sd1;
    out += alpha_i * alpha_i * stats::dnorm(m1, m1, sqrt(2) * sd1);
    for (int j=i+1; j < N1; j++) {
      out += 2 * alpha_i * weights1(j) * stats::dnorm(
        m1, params1(j, 0), sqrt(var1 + std::pow(params1(j, 1), 2)));
    }
  }

  for (int i=0; i < N2; i++) {
    double alpha_i = weights2(i);
    double m1 = params2(i, 0);
    double sd1 = params2(i, 1);
    double var1 = sd1 * sd1;
    out += alpha_i * alpha_i * stats::dnorm(m1, m1, sqrt(2) * sd1);
    for (int j=i+1; j < N1; j++) {
      out += 2 * alpha_i * weights2(j) * stats::dnorm(
        m1, params2(j, 0), sqrt(var1 + std::pow(params2(j, 1), 2)));
    }
  }

  for (int i=0; i < N1; i++) {
    double alpha_i = weights1(i);
    double m1 = params1(i, 0);
    double sd1 = params1(i, 1);
    double var1 = sd1 * sd1;
    out += alpha_i * alpha_i * stats::dnorm(m1, m1, sqrt(2) * sd1);
    for (int j=0; j < N2; j++) {
      out -= 2 * alpha_i * weights2(j) * stats::dnorm(
        m1, params2(j, 0), sqrt(var1 + std::pow(params2(j, 1), 2)));
    }
  }
  return out;
}

template<typename Hierarchy>
void MetropolisedMarginalSemiHdp<Hierarchy>::computeDistances() {
  for (int g=0; g < numGroups; g++) {
    std::tuple<arma::mat, arma::vec> paramsAndWeights1 = getParamsAndWeights(g);
    #pragma omp parallel for
    for (int k=g+1; k < numGroups; k++) {
      std::tuple<arma::mat, arma::vec> paramsAndWeights2 = getParamsAndWeights(k);
      distances(g, k) = gm_distance(
        std::get<0>(paramsAndWeights1), std::get<1>(paramsAndWeights1),
        std::get<0>(paramsAndWeights2), std::get<1>(paramsAndWeights2));
      distances(k, g) = distances(g, k);
    }
  }
}


template<typename Hierarchy>
void MetropolisedMarginalSemiHdp<Hierarchy>::sampleGroupAssignments() {
  computeDistances();
  for (int g=0; g < numGroups; g++) {
    int prev = c(g);
    arma::vec dists = distances.row(g) + 1.0;
    arma::vec probas = arma::pow(dists, -1);
    probas /= (2.0 *arma::max(probas));
    probas += 0.5;
    int proposed = stats::rdiscrete(probas);
    double logProbaCurr = std::log(omegas(c(g)));
    double logProbaProposed = std::log(omegas(c(proposed)));
    for (int j=0; j < samplesPerGroup(g); j++) {
      logProbaCurr += logLikeForComponent(data[g][j], c(g));
      logProbaProposed += logLikeForComponent(data[g][j], proposed);
    }
    // TODO add prior over C!
    double logratio = logProbaProposed - logProbaCurr;
    double arate = std::min(1.0, exp(logratio));
    if (stats::runif(0.0, 1.0, engine) < arate)
      c(g) = proposed;

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


#endif  // CPPMODEL_METROPOLISED_MARGINAL_SEMI_HDP_IMP_HPP
