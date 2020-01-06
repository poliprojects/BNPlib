//
// Created by mario on 20/05/19.
//

#ifndef CPPMODEL_METROPOLIS_HASTINGS_HPP
#define CPPMODEL_METROPOLIS_HASTINGS_HPP

#include "options.hpp"
#include "random_engine.hpp"


/*
 * Abstract class for Metropolis Hastings sampling
 */

template<typename Theta>
class MetropolisHastings {
 protected:
  double acceptanceRate;
  int numSteps;
  Theta currState;
  stats::rand_engine_t& engine = RandomEngine::Instance().get();

 public:
  virtual ~MetropolisHastings() {};

  MetropolisHastings() {}

  MetropolisHastings(Theta state): currState(state) {
    acceptanceRate = 0.0;
    numSteps = 0;
  }

  virtual Theta propose() = 0;

  virtual double logratio(Theta newState) = 0;

  void step() {
    Theta newState = propose();
    double r = stats::runif(0.0, 1.0, engine);
    double arate = std::min(1.0, exp(logratio(newState)));
    numSteps += 1;
    acceptanceRate += (arate - acceptanceRate) / numSteps;
    if (r < arate) {
      currState = newState;
    }
  }

  void resetState(Theta state) {
    currState = state;
    numSteps = 0;
    acceptanceRate = 0.0;
  }

  double getAcceptanceRate() const {
    return acceptanceRate;
  }

  const Theta &getState() {
    return currState;
  }

  void setCurrState(Theta currState) {
    MetropolisHastings::currState = currState;
  }

};

#endif //CPPMODEL_METROPOLIS_HASTINGS_HPP
