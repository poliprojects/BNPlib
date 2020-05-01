#ifndef HIERARCHYNNIG_IMP_HPP
#define HIERARCHYNNIG_IMP_HPP

#include "HierarchyNNIG.hpp"


template<class Hypers> 
void HierarchyNNIG<Hypers>::check_state_validity(){
    assert(this->state[1](0,0) > 0);
}


template<class Hypers> 
Eigen::VectorXd HierarchyNNIG<Hypers>::like(const Eigen::MatrixXd &data){
    Eigen::VectorXd result(data.rows());
    for(size_t i = 0; i < data.rows(); i++){
        result(i) = exp(stan::math::normal_lpdf(data(i,0), this->state[0](0,0),
            this->state[1](0,0)));
    }
    return result;
}


template<class Hypers> 
void HierarchyNNIG<Hypers>::draw(){
    double sigma2_new = sqrt( stan::math::inv_gamma_rng(
        this->hypers->get_alpha0(), this->hypers->get_beta0(), this->rng) );

    double mu_new = stan::math::normal_rng(this->hypers->get_mu0(),
        sqrt(sigma2_new/this->hypers->get_lambda()), this->rng);
    
    this->state[0](0,0) = mu_new;
    this->state[1](0,0) = sigma2_new;
}


template<class Hypers> 
Eigen::VectorXd HierarchyNNIG<Hypers>::eval_marg(const Eigen::MatrixXd &data){

    double sigtilde = sqrt( this->hypers->get_beta0()*(
        this->hypers->get_lambda()+1) / ( this->hypers->get_alpha0() *
        this->hypers->get_lambda() ) );
   
    Eigen::VectorXd result(data.rows());
    for(size_t i = 0; i < data.rows(); i++){
        result(i) = exp( stan::math::student_t_lpdf(data(i,0),
            2*this->hypers->get_alpha0(), this->hypers->get_mu0(), sigtilde) );
    }
    return result;
}


template<class Hypers> 
std::vector<double> HierarchyNNIG<Hypers>::normal_gamma_update(
    const Eigen::VectorXd &data, const double mu0, const double alpha0,
    const double beta0, const double lambda0){

    double mu_post, alpha_post, beta_post, lambda_post;
    unsigned int n = data.rows();

    if(n == 0){
        return std::vector<double>{mu0, alpha0, beta0, lambda0};
    }
    
    // Compute sample mean
    double y_bar = data.mean();

    // Compute parameters
    mu_post = (lambda0 * mu0 + n * y_bar) / (lambda0 + n);
    alpha_post = alpha0 + 0.5 * n;
    // double ss = n * arma::var(data, 1); // divides by n, not n-1
    double ss = (data.dot(data)) - n*y_bar*y_bar;
    beta_post = beta0 + 0.5 * ss + 0.5 * lambda0 * n * std::pow((
        y_bar - mu0), 2) / (n + lambda0);
    lambda_post = lambda0 + n;
    
    return std::vector<double>{mu_post, alpha_post, beta_post, lambda_post};
}


template<class Hypers> 
void HierarchyNNIG<Hypers>::sample_given_data(const Eigen::MatrixXd &data){

    // Get current values of parameters
    double mu0     = this->hypers->get_mu0();
    double lambda0 = this->hypers->get_lambda();
    double alpha0  = this->hypers->get_alpha0();
    double beta0   = this->hypers->get_beta0();

    std::vector<double> temp = normal_gamma_update(data.col(0), mu0, alpha0,
        beta0, lambda0);

    double mu_post     = temp[0];
    double alpha_post  = temp[1];
    double beta_post   = temp[2];
    double lambda_post = temp[3];

    // Get a sample
    double sigma2_new,mu_new;
    sigma2_new = sqrt(stan::math::inv_gamma_rng(alpha_post, beta_post,
        this->rng));
    mu_new = stan::math::normal_rng(mu_post, sqrt(sigma2_new/lambda_post),
        this->rng); 
    this->state[0](0,0) = mu_new;
    this->state[1](0,0) = sigma2_new;
}


#endif // HIERARCHYNNIG_IMP_HPP
