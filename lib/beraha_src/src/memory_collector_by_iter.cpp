#include "memory_collector_by_iter.hpp"


void MemoryCollectorByIter::collect(std::string paramName, arma::mat val) {
  if (currStep == 0)
    chains.resize(numSteps);
  chains[currStep][paramName] = val;
}

void MemoryCollectorByIter::setMaxNumSteps(int maxNumSteps) {
  MemoryCollectorByIter::maxNumSteps = maxNumSteps;
}

void MemoryCollectorByIter::setNumSteps(int numSteps) {
  MemoryCollectorByIter::numSteps = numSteps;
}
