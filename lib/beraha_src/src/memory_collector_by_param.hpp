//
// Created by mario on 22/05/19.
//

#ifndef CPPMODEL_MEMORY_COLLECTOR_BYPARAM_HPP
#define CPPMODEL_MEMORY_COLLECTOR_BYPARAM_HPP

#include "options.hpp"
#include <sstream>
#include <fstream>

#ifdef USE_R
#include <Rcpp.h>
#endif

/*
 * Collects in memory *ALL* the parameters of the MCMC simulation
 * by storing a single chain for each parameter. Useful when the
 * state size (number of parameters) is fixed
 */
class MemoryCollectorByParam {
 protected:
  std::map<std::string, std::vector<arma::mat>> chains;
  int maxNumSteps;
  int numSteps;


 public:
  ~MemoryCollectorByParam() {}

  MemoryCollectorByParam() {}
  MemoryCollectorByParam(int numSteps): maxNumSteps(numSteps), numSteps(numSteps) {}

  void collect(std::string paramName, arma::mat val);

  void collect(std::string paramName, double val);

  void to_csv(std::string fileName);

  int getNumSteps() const;

  arma::mat get(std::string paramName, int iternum) const;

  arma::mat getParamChain(std::string paramName) {
    arma::mat out(chains[paramName].size(), chains[paramName][0].n_elem);
    for (int i=0; i < chains[paramName].size(); i++) {
      out.row(i) = chains[paramName][i].t();
    }
    return out;
  }

  void setMaxNumSteps(int maxNumSteps);
  void setNumSteps(int numSteps);

  #ifdef USE_R
    Rcpp::List getChains() {
      Rcpp::List out = Rcpp::List::create();
      for (auto it: chains) {
        arma::field<arma::mat> curr(it.second.size());
        for (int i = 0; i < it.second.size(); i++)
          curr(i) = it.second[i];

        out[it.first] = Rcpp::wrap(curr);
      }
      return out;
    }

    void restoreFromR(const Rcpp::List &Rchains) {
      std::vector<std::string> names = Rchains.names();
      for (std::string name: names) {
        chains[name] = Rcpp::as<std::vector<arma::mat>>(Rchains[name]);
      }
      numSteps = chains[names[0]].size();
      maxNumSteps = numSteps;
    }
  #endif
};

#endif //CPPMODEL_MEMORY_COLLECTOR_BYPARAM_HPP
