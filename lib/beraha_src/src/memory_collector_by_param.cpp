//
// Created by mario on 22/05/19.
//

#include "memory_collector_by_param.hpp"

void MemoryCollectorByParam::collect(std::string paramName, arma::mat val) {
  if (chains.find(paramName) == chains.end())
    chains[paramName].reserve(maxNumSteps + 10);

  chains[paramName].push_back(val);
}

void MemoryCollectorByParam::collect(std::string paramName, double val) {
  collect(paramName, arma::vec{val});
}

void MemoryCollectorByParam::to_csv(std::string fileName) {
  numSteps = chains.begin()->second.size();

  std::ostringstream firstRow;
  for (auto it: chains) {
    std::string name = it.first;
    int n_rows = it.second[0].n_rows;
    int n_cols = it.second[0].n_cols;

    if ((n_rows == 1) || (n_cols == 1)) {
      for (int i = 0; i < it.second[0].n_elem; i++)
        firstRow << name << "." << i + 1 << ", ";
    } else {

      for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
          firstRow << name << "." << i + 1 << "." << j + 1 << ", ";
        }
      }
    }
  }
  std::string firstRowStr = firstRow.str();

  // firstRow ends with ", " --> remove last 2 chars
  firstRowStr.pop_back();
  firstRowStr.pop_back();

  // erase file if already present
  std::ofstream ofs;
  ofs.open(fileName, std::ofstream::out | std::ofstream::trunc);
  ofs.close();

  std::ofstream fp(fileName);
  if (fp.is_open()) {
    fp << firstRowStr << "\n";

    for (int i = 0; i < numSteps; i++) {
      arma::rowvec currRow;
      for (auto it: chains)
        currRow = arma::join_rows(
            currRow, arma::vectorise(it.second[i], 1));

      currRow.save(fp, arma::csv_ascii);
    }
  }
}

arma::mat MemoryCollectorByParam::get(std::string paramName, int iternum) const {
  assert(iternum <= numSteps);
  return chains.at(paramName)[iternum];
}


int MemoryCollectorByParam::getNumSteps() const {
  return numSteps;
}
void MemoryCollectorByParam::setMaxNumSteps(int maxNumSteps) {
  MemoryCollectorByParam::maxNumSteps = maxNumSteps;
}
void MemoryCollectorByParam::setNumSteps(int numSteps) {
  MemoryCollectorByParam::numSteps = numSteps;
}
