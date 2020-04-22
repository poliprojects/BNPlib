#ifndef HIERARCHYNW_IMP_HPP
#define HIERARCHYNW_IMP_HPP

#include "HierarchyNW.hpp"


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
    this->state[1] = stan::math::multi_normal_rng( this->hypers->get_nu(),
    	???, this->rng )
    //double sigma2_new = stan::math::inv_gamma_rng(this->hypers->get_alpha0(),
    //    this->hypers->get_beta0(), this->rng);
}


template<class Hypers> 
Eigen::VectorXd HierarchyNW<Hypers>::eval_marg(const Eigen::MatrixXd &data){
    Eigen::VectorXd result(data.cols());
    for(int i = 0; i < data.cols(); i++){
        result(i) = data(0,i) / 20.0;
    }
    return result;
}


template<class Hypers> 
std::vector<Eigen::MatrixXd> HierarchyNW<Hypers>::dummy_update(
    const Eigen::MatrixXd &data, const Eigen::VectorXd &mu0, const
    Eigen::MatrixXd &lambda0){

    int n = data.cols();
    Eigen::VectorXd mu_post = mu0 * std::pow(0.99,n);
    Eigen::MatrixXd lambda_post = lambda0 * std::pow(0.99,n);
    
    return std::vector<Eigen::MatrixXd>{mu_post, lambda_post};
}


template<class Hypers> 
void HierarchyNW<Hypers>::sample_given_data(const Eigen::MatrixXd &data){
    // Get current values of parameters
    Eigen::VectorXd mu0 = this->hypers->get_mu0();
    Eigen::MatrixXd lambda0 = this->hypers->get_lambda0();

    std::vector<Eigen::MatrixXd> temp = dummy_update(data, mu0, lambda0);

    Eigen::VectorXd mu_post = temp[0];
    Eigen::MatrixXd lambda_post = temp[1];

    // Get a sample
    this->state[0] = mu_post * 0.5;
    this->state[1] = lambda_post * 0.5;
}


#endif // HIERARCHYNW_IMP_HPP
