//
// Created by mario on 16/05/19.
//

#ifndef CPPMODEL_UTILS_HPP
#define CPPMODEL_UTILS_HPP

#include "options.hpp"

double logSumExp(arma::vec logprobas);

// Returns a vector of postMean, postA, postB, postLambda
arma::vec normalGammaUpdate(
    arma::vec data, double priorMean, double priorA, double priorB,
    double priorLambda);

template<typename T>
arma::vec stdToArma(std::vector<T> vec) {
  arma::vec out(vec.size());
  for (int i = 0; i < vec.size(); i++) {
    out(i) = (double) vec[i];
  }
  return out;
}

template<class Matrix>
void print_matrix(Matrix matrix) {
  matrix.print(std::cout);
}

arma::vec twoNormalMixture(int numSamples, double mean1, double sd1,
                           double mean2, double sd2, double w);


arma::vec uniformNormalMixture(int numSamples, arma::vec means, arma::vec sds) ;

//provide explicit instantiations of the template function for
//every matrix type you use somewhere in your program.
template void print_matrix<arma::mat>(arma::mat matrix);
template void print_matrix<arma::cx_mat>(arma::cx_mat matrix);

#endif //CPPMODEL_UTILS_HPP
