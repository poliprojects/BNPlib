#ifndef CPPMODEL_MEMORY_COLLECTOR_BYITER_HPP
#define CPPMODEL_MEMORY_COLLECTOR_BYITER_HPP

#include "options.hpp"
#include <sstream>
#include <fstream>

#ifdef USE_R
#include <Rcpp.h>
#endif

/*
 * Collects in memory *ALL* the parameters of the MCMC simulation
 */
class MemoryCollectorByIter {
 protected:
  std::vector<std::map<std::string, arma::mat>> chains;
  int maxNumSteps;
  int numSteps;
  int currStep = -1;


 public:
  ~MemoryCollectorByIter() {}

  MemoryCollectorByIter() {}
  
  MemoryCollectorByIter(int numSteps):
      maxNumSteps(numSteps), numSteps(maxNumSteps) {
    chains.resize(numSteps);
  }

  void collect(std::string paramName, arma::mat val);

  void newIter() {currStep += 1;}

  arma::mat get(std::string paramName, int iternum) const {
    return chains.at(iternum).at(paramName);
  }

  arma::mat getParamChain(std::string paramName) {
    // !!! Just for debugging, might crash if the parameter size
    // chanegs across iterations !!!
    std::cout << "Getting param: " << paramName << std::endl;
    arma::mat out(chains.size(), chains[0][paramName].n_elem);
    std::cout << "Initialized matrix" << std::endl;
    for (int i=0; i < chains.size(); i++) {
      out.row(i) = chains[i][paramName].t();
    }
    return out;
  }
  
  int getNumSteps() const {
    return numSteps;
  }

  void setMaxNumSteps(int maxNumSteps);
  void setNumSteps(int numSteps);

  #ifdef USE_R
    Rcpp::List getChains() {
      Rcpp::List out = Rcpp::List(currStep);
      for (int i = 0; i < currStep; i++) {
        Rcpp::List iterList = Rcpp::List::create();
        for (auto it: chains[i]) {
          iterList[it.first] = Rcpp::wrap(it.second);
        }
        out[i] = iterList;
      }
      return out;
    }

    void restoreFromR(const Rcpp::List &Rchains) {
      numSteps = Rchains.size();
      maxNumSteps = numSteps;
      chains.resize(numSteps);
      for (int i = 0; i < numSteps; i++) {
        std::map<std::string, arma::mat> iter;
        Rcpp::List currIter = Rchains[i];
        std::vector<std::string> names = currIter.names();
        for (std::string name: names)
          iter[name] = Rcpp::as<arma::mat>(currIter[name]);
        chains[i] = iter;
      }
    }
  #endif
};

#endif  // CPPMODEL_MEMORY_COLLECTOR_BYITER_HPP
