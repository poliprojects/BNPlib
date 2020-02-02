#ifndef NNIGHIERARCHY_HPP
#define NNIGHIERARCHY_HPP

#include <array>
#include <random>
#include <vector>
#include <memory>
#include <stan/math/prim/mat.hpp>

#include "includes_universal.hpp"

template<class Hypers> //Hypers = TupleWrapper, distro, ...
class NNIGHierarchy {
protected:
    using state_t = std::array<par_t,2>;

    std::mt19937 rng;
    state_t state; // current values for F's parameters: mu, sigma
    std::shared_ptr<Hypers> hypers; // current values for G0's parameters:
                                    // mu_0,Lambda0, alpha, beta

public:
    // Contructors:
    ~NNIGHierarchy() = default;
    NNIGHierarchy(std::shared_ptr<Hypers> hypers): hypers(hypers) {}

    // Getters/setters:
    state_t get_state(){return state;}
    void set_state(const state_t &s){state = s;}
    void set_state(int pos, par_t val){state[pos] = val;}

    int get_count(){return hypers.use_count();}

    double log_like(data_t datum){
        return exp(stan::math::normal_lpdf(datum, state[0], state[1]));
    }

    void draw(){
    float sigma_new = stan::math::inv_gamma_rng(hypers->get_alpha0(),
        hypers->get_beta0(), rng);
    float mu_new = stan::math::normal_rng(hypers->get_mu0(),
        sigma_new/hypers->get_lambda(), rng);
    state[0] = mu_new;
    state[1] = sigma_new;
    }

    void sample_given_data(std::vector<data_t> data){
    // Get current values of parameters
    auto mu0     = hypers->get_mu0();
    auto lambda0 = hypers->get_lambda();
    auto alpha0  = hypers->get_alpha0();
    auto beta0   = hypers->get_beta0();
    parvec_t temp = normal_gamma_update(data, mu0, alpha0, beta0, lambda0);
    auto mu_post = temp[0];
    auto alpha_post = temp[1];
    auto beta_post = temp[2];
    auto lambda_post = temp[3];
    // Get a sample
    par_t sigma_new = stan::math::inv_gamma_rng(alpha_post, beta_post, rng);
    par_t mu_new = stan::math::normal_rng(mu_post, sigma_new/lambda_post,
        rng); //? is it ok /lambda_post?
    state[0] = mu_new;
    state[1] = sigma_new;
  }


  parvec_t normal_gamma_update(std::vector<data_t> data,
    double mu0, double alpha0, double beta0, double lambda0) {
    
    double mu_post, alpha_post, beta_post, lambda_post;
    int n = data.size();
    if (n == 0)
        return parvec_t{mu0, alpha0, beta0, lambda0};
    
    // Compute sample mean
    double y_bar = accumulate(data.begin(), data.end(), 0.0) / n;
    mu_post = (lambda0 * mu0 + n * y_bar) / (lambda0 + n);
    alpha_post = 1.0 * alpha0 + 1.0 * n / 2;
    // double ss = n * arma::var(data, 1); // divides by n, not n-1
    double ss = std::inner_product( data.begin(), data.end(), data.begin(),
        0.0 ) - n*y_bar*y_bar; // TODO check
    beta_post = (beta0 + 0.5 * ss + 0.5 * lambda0 /
        (n + lambda0) * n * std::pow((y_bar - mu0), 2));
    lambda_post = lambda0 + n;
    return parvec_t{mu_post, alpha_post, beta_post, lambda_post};
}

};

#endif // NNIGHIERARCHY_HPP
