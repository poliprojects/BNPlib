#ifndef HIERARCHYNW_IMP_HPP
#define HIERARCHYNW_IMP_HPP

#include "HierarchyNW.hpp"

// TODO check everything

template<class Hypers> 
Eigen::VectorXd HierarchyNW<Hypers>::like(const Eigen::MatrixXd &data){
    Eigen::VectorXd result(data.rows());
    for(int i = 0; i < data.rows(); i++){
        result(i) = exp( stan::math::multi_normal_lpdf(data(i), this->state[0],
            inverse) );
    }
    return result;
}


template<class Hypers> 
void HierarchyNW<Hypers>::draw(){
    this->state[0] = stan::math::multi_normal_rng( this->hypers->get_mu0(),
        inverse*(1/this->hypers->get_lambda()), this->rng );
    Eigen::MatrixXd sigma_new = stan::math::wishart_rng( this->hypers->get_nu(),
    	this->hypers->get_tau0_inv(), this->rng );
    this->state[1] = sigma_new.inverse();
}


template<class Hypers> 
Eigen::VectorXd HierarchyNW<Hypers>::eval_marg(const Eigen::MatrixXd &data){
    Eigen::VectorXd result(data.cols());
    // TODO da fare
    return result;
}


template<class Hypers> 
std::vector<Eigen::MatrixXd> HierarchyNW<Hypers>::normal_wishart_update(
    const Eigen::MatrixXd &data, const Eigen::VectorXd &mu0,
    const double lambda, const Eigen::MatrixXd &tau0, const double nu){
    int n = data.cols();
    Eigen::MatrixXd lambda_post(1,1), nu_post(1,1);
    Eigen::VectorXd mubar = data.rowwise().mean(); // TODO ??? sparato a caso

    lambda_post(0,0) = lambda + n;
    nu_post(0,0) = nu + n;
    Eigen::VectorXd mu_post = (lambda*mu0 + n*mubar) * (1/lambda+n);

    // Compute tau_post
    Eigen::MatrixXd tau_post = tau0;
    for(unsigned int i = 0; i < n; i++){
    	tau_post = tau_post + (data(i)-mubar)*(data(i)-mubar).transpose();
    }
    tau_post = tau_post + (nu*lambda/(nu+lambda)) *
    	(mubar-mu0)*(mubar-mu0).transpose();
    
    return std::vector<Eigen::MatrixXd>{mu_post,lambda_post,tau_post,nu_post};
}


template<class Hypers> 
void HierarchyNW<Hypers>::sample_given_data(const Eigen::MatrixXd &data){
    // Get current values of parameters
    Eigen::VectorXd mu0 = this->hypers->get_mu0();
    double lambda = this->hypers->get_lambda();
    Eigen::MatrixXd tau0 = this->hypers->get_tau0();
    double nu = this->hypers->get_nu();

    std::vector<Eigen::MatrixXd> temp = normal_wishart_update(data, mu0, lambda,
    	tau0, nu);

    Eigen::VectorXd mu_post = temp[0];
    double lambda_post = temp[1];
    Eigen::MatrixXd tau_post = temp[2];
    double nu = temp[3];

    sigma_post = tau_post.inverse();

    // Get a sample
    Eigen::VectorXd mu_new = stan::math::multi_normal_rng(mu_post,
        sigma_post*(1/lambda_post), this->rng);
    Eigen::MatrixXd sigma_new = stan::math::wishart_rng(nu_post, sigma_post,
    	this->rng);
    this->state[0] = mu_new;
    this->state[1] = sigma_new.inverse();
}


#endif // HIERARCHYNW_IMP_HPP
