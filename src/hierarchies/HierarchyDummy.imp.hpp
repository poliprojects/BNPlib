#ifndef HIERARCHYDUMMY_IMP_HPP
#define HIERARCHYDUMMY_IMP_HPP

#include "HierarchyDummy.hpp"



template<class Hypers> 
Eigen::VectorXd HierarchyDummy<Hypers>::like(const Eigen::MatrixXd &datum){
    Eigen::VectorXd result(datum.cols());
    for(int i = 0; i < datum.cols(); i++){
        result(i) = datum(0,i) / 500.0;
    }
    return result;
}

template<class Hypers> 
void HierarchyDummy<Hypers>::draw(){
    for(int i = 0; i < state[0].size(); i++){
        state[0](i) = 2;
    }
    state[1].setIdentity();
    //double sigma2_new = stan::math::inv_gamma_rng(hypers->get_alpha0(),
    //    hypers->get_beta0(), rng);
}


template<class Hypers> 
Eigen::VectorXd HierarchyDummy<Hypers>::eval_marg(const Eigen::MatrixXd &datum){
    Eigen::VectorXd result(datum.cols());
    for(int i = 0; i < datum.cols(); i++){
        result(i) = datum(0,i) / 1000.0;
    }
    return result;
}


template<class Hypers> 
std::vector<Eigen::MatrixXd> HierarchyDummy<Hypers>::dummy_update(
    const Eigen::MatrixXd &data, const Eigen::VectorXd &mu0, const
    Eigen::MatrixXd &lambda0){

    int n = data.cols();
    Eigen::VectorXd mu_post = mu0 * std::pow(0.99,n);
    Eigen::MatrixXd lambda_post = lambda0 * std::pow(0.99,n);
    
    return std::vector<Eigen::MatrixXd>{mu_post, lambda_post};
}


template<class Hypers> 
void HierarchyDummy<Hypers>::sample_given_data(const Eigen::MatrixXd &data){
    // Get current values of parameters
    Eigen::VectorXd mu0 = hypers->get_mu0();
    Eigen::MatrixXd lambda0 = hypers->get_lambda0();

    std::vector<Eigen::MatrixXd> temp = dummy_update(data, mu0, lambda0);

    Eigen::VectorXd mu_post = temp[0];
    Eigen::MatrixXd lambda_post = temp[1];

    // Get a sample
    state[0] = mu_post * 0.5;
    state[1] = lambda_post * 0.5;
}


#endif // HIERARCHYDUMMY_IMP_HPP
