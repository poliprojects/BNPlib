//
// Created by mario on 10/05/19.
//

#ifndef CPPMODEL_CLUSTER_HPP
#define CPPMODEL_CLUSTER_HPP

#include "options.hpp"

template<typename Theta>
class Cluster {
 protected:
  Theta param;
  vector<double> data;
  arma::uvec nDataFromRestaurant;
  arma::uvec nPseudoDataFromRestaurant;
  int numElems = 0;
  int numPseudoElems = 0;

 public:
  ~Cluster() = default;

  Cluster<Theta>() {data.reserve(10);}

  Cluster<Theta>(Theta param): param(param) {
    data.reserve(10);
    // We assume there will be at most 100 groups
    nDataFromRestaurant = arma::uvec(100, arma::fill::zeros);
    nPseudoDataFromRestaurant = arma::uvec(100, arma::fill::zeros);
    param.samplePrior();
  }

  Cluster<Theta>(Theta param, int ngroups): param(param) {
    data.reserve(10);
    nDataFromRestaurant = arma::uvec(ngroups, arma::fill::zeros);
    nPseudoDataFromRestaurant = arma::uvec(ngroups, arma::fill::zeros);
    param.samplePrior();
  }

  Cluster<Theta>(int ngroups) {
    data.reserve(10);
    nDataFromRestaurant = arma::uvec(ngroups, arma::fill::zeros);
    nPseudoDataFromRestaurant = arma::uvec(ngroups, arma::fill::zeros);
    param.samplePrior();
  }

  void add(double datum) {
    data.push_back(datum);
    numElems += 1;
  }

  void add(double datum, int res, bool resample=false) {
    data.push_back(datum);
    nDataFromRestaurant(res) += 1;
    numElems += 1;
    if (resample)
      sample();
  }

  bool isEmpty() { return data.empty(); }

  Theta getParam() const {
    return param;
  }

  Theta getParamForMerge() {
    return param;
  }

  const vector<double> &getData() const {
    return data;
  }

  void sample() {
    if (data.empty())
      param.samplePrior();
    else
      param.sample(data);
  }

  double loglike(double x) {
    return param.loglike(x);
  }

  void remove(double datum) {
    auto posIt = std::find(data.begin(), data.end(), datum);
    if (posIt != data.end()) {
      std::swap(*posIt, data.back());
      data.pop_back();
      numElems -= 1;
    } else {
      std::cout << "Tried to remove datum " << datum << " but did not find it" << std::endl;
      std::cout << "Current data: " << std::endl;
      for (double d : data) {
        std::cout << d << ", ";
      }
      std::cout << std::endl;
    }
    if (data.empty())
      data.reserve(10);
    assert(numElems >= 0);
  }

  void remove(double datum, int res) {
    remove(datum);
    nDataFromRestaurant(res) -= 1;
    assert(nDataFromRestaurant(res)>=0);
  }

  void empty() {
    data.clear();
    param.reInitialize();
  }

  void merge(Cluster<Theta>& other) {
    if (! other.isEmpty()) {
      std::vector<double> otherData = other.getData();
      data.reserve(data.size() + otherData.size() + 10);
      data.insert(
          data.end(), otherData.begin(), otherData.end());
      numElems = data.size();
      nDataFromRestaurant = other.getNumDataFromRestaurant();
      setParam(other.getParamForMerge());
    }
  }

  void setParam(Theta param) {
    Cluster::param = param;
  }

  int size(bool usePseudo=false) {
    int out;
    if (!usePseudo) {
      out = data.size() == 0 ? numElems : data.size();
    }
    else
      out = numElems == 0 ? numPseudoElems : numElems;
    return out;
  }

  int size(int res) {
    int out;
    if (arma::sum(nDataFromRestaurant) == 0)
      out = nPseudoDataFromRestaurant(res);
    else
      out = nDataFromRestaurant(res);
    return out;
  }

  /*
   * Just for Neal's algorithm 2
   */
  double logProbaForDatum(double x) {
    int n = numElems == 0 ? numPseudoElems : numElems;
    double out = log(1.0 * n) + param.loglike(x);
    return out;
  }

  arma::uvec getNumDataFromRestaurant() {
    return nDataFromRestaurant;
  }

  void setNumElems(int numElems) {
    Cluster::numElems = numElems;
  }

  void setNumPseudoElems(int numElems) {
    Cluster::numPseudoElems = numElems;
  }

  void setNDataFromRestaurant(arma::uvec &nDataFromRestaurant) {
    Cluster::nDataFromRestaurant = nDataFromRestaurant;
    Cluster::numElems = arma::sum(nDataFromRestaurant);
  }

  void setNPseudoDataFromRestaurant(arma::uvec &nDataFromRestaurant) {
    Cluster::nPseudoDataFromRestaurant = nDataFromRestaurant;
    Cluster::numPseudoElems = arma::sum(nDataFromRestaurant);
  }

};

#endif //CPPMODEL_CLUSTER_HPP
