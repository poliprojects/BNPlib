#ifndef CPPMODEL_RCPP_FUNCTIONS_H
#define CPPMODEL_RCPP_FUNCTIONS_H

// [[Rcpp::depends(RcppArmadillo)]]
#include "options.hpp"
#include "blocked_gibbs.hpp"
#include "semi_hdp_blocked_gibbs.hpp"
#include "marginal_semi_hdp.hpp"
#include "memory_collector_by_param.hpp"
#include "memory_collector_by_iter.hpp"

// [[Rcpp::export]]
Rcpp::List runGibbs(std::string model, const Rcpp::List& data_, int numIter,
                    int burnIn, int warmup, int thin, int numComponents=-1,
                    int log_every=100, bool usePseudoPrior=false);

// [[Rcpp::export]]
arma::mat predictiveSamples(std::string model, const Rcpp::List &chains);

// [[Rcpp::export]]
Rcpp::List runNormalConjugateGibbs(
    const Rcpp::List& data_, int numIter, int burnIn, int warmup, int thin,
    int numComponents, int log_every, bool usePseudoPrior);

// [[Rcpp::export]]
Rcpp::List runSemiHdpBlockedGibbs(
    const Rcpp::List& data_, int numIter, int burnIn, int warmup, int thin,
    int numComponents, int log_every, bool usePseudoPrior);

// [[Rcpp::export]]
Rcpp::List runSemiHdpMarginal(
    const Rcpp::List& data_, int numIter, int burnIn, int warmup, int thin,
    int log_every, bool usePseudoPrior);




#endif //CPPMODEL_RCPP_FUNCTIONS_H
