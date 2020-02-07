#ifndef NNIGHIERARCHY_IMP_HPP
#define NNIGHIERARCHY_IMP_HPP

#include "NNIGHierarchy.hpp"

template<class Hypers> 
double NNIGHierarchy<Hypers>::log_like(double datum){
        return exp(stan::math::normal_lpdf(datum, state[0], state[1]));
}

template<class Hypers> 
Eigen::VectorXf NNIGHierarchy<Hypers>::log_like(std::vector<double> datum){
    // TODO stan per vector?? 
    Eigen::VectorXf result;
    for (int i=0; i< datum.size(); i++)
    result(i) = exp(stan::math::normal_lpdf(datum[i], state[0], state[1]));
    
    return result;
}

template<class Hypers> 
void NNIGHierarchy<Hypers>::draw(){
    float sigma_new = stan::math::inv_gamma_rng(hypers->get_alpha0(),
        hypers->get_beta0(), rng);
    float mu_new = stan::math::normal_rng(hypers->get_mu0(),
        sigma_new/hypers->get_lambda(), rng);
    state[0] = mu_new;
    state[1] = sigma_new;
    }

template<class Hypers> 
Eigen::VectorXf NNIGHierarchy<Hypers>::eval_G0(std::vector<double> datum){
     // TODO stan per vector??
    Eigen::VectorXf result;
    for (int i=0; i< datum.size(); i++)
    result(i)=exp(stan::math::normal_lpdf(datum[i] , hypers->get_mu0(),
        stan::math::inv_gamma_lpdf(datum[i],hypers->get_alpha0(),
        hypers->get_beta0())/hypers->get_lambda() ) );
    
    return result;
    }

template<class Hypers> 
double NNIGHierarchy<Hypers>::eval_G0(double datum){ 
    return exp(stan::math::normal_lpdf(datum , hypers->get_mu0(),
        stan::math::inv_gamma_lpdf(datum,hypers->get_alpha0(),
        hypers->get_beta0())/hypers->get_lambda() ) ); 
    }

template<class Hypers> 
void NNIGHierarchy<Hypers>::sample_given_data(std::vector<double> data){
    // Get current values of parameters
    auto mu0     = hypers->get_mu0();
    auto lambda0 = hypers->get_lambda();
    auto alpha0  = hypers->get_alpha0();
    auto beta0   = hypers->get_beta0();

    std::vector<double> temp = normal_gamma_update(data, mu0, alpha0, beta0,
        lambda0);

    auto mu_post = temp[0];
    auto alpha_post = temp[1];
    auto beta_post = temp[2];
    auto lambda_post = temp[3];

    // Get a sample
    double sigma_new = stan::math::inv_gamma_rng(alpha_post, beta_post, rng);
    double mu_new = stan::math::normal_rng(mu_post, sigma_new/lambda_post,
        rng); //? is it ok /lambda_post?
    state[0] = mu_new;
    state[1] = sigma_new;
  }


template<class Hypers> 
std::vector<double> NNIGHierarchy<Hypers>::normal_gamma_update(
    std::vector<double> data, double mu0, double alpha0, double beta0,
    double lambda0) {

    double mu_post, alpha_post, beta_post, lambda_post;
    int n = data.size();

    if (n == 0){
        return std::vector<double>{mu0, alpha0, beta0, lambda0};
    }
    
    // Compute sample mean
    double y_bar = accumulate(data.begin(), data.end(), 0.0) / n;

    // Compute parameters
    mu_post = (lambda0 * mu0 + n * y_bar) / (lambda0 + n);
    alpha_post = 1.0 * alpha0 + 1.0 * n / 2;
    // double ss = n * arma::var(data, 1); // divides by n, not n-1
    double ss = std::inner_product( data.begin(), data.end(), data.begin(),
        0.0 ) - n*y_bar*y_bar; // TODO check
    beta_post = (beta0 + 0.5 * ss + 0.5 * lambda0 /
        (n + lambda0) * n * std::pow((y_bar - mu0), 2));
    lambda_post = lambda0 + n;
    
    return std::vector<double>{mu_post, alpha_post, beta_post, lambda_post};
}

#endif // NNIGHIERARCHY_IMP_HPP
