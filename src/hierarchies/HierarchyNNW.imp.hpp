#ifndef HIERARCHYNNW_IMP_HPP
#define HIERARCHYNNW_IMP_HPP

#include "HierarchyNNW.hpp"


template<class Hypers> 
Eigen::VectorXd HierarchyNNW<Hypers>::like(const Eigen::MatrixXd &data){
    unsigned int n = data.rows();
    Eigen::VectorXd result(n);
	Eigen::MatrixXd sigma = this->state[1].inverse();
	EigenRowVec mu(this->state[0]);
    for(unsigned int i = 0; i < n; i++){
        EigenRowVec datum = data.row(i);
        result(i) = std::exp(stan::math::multi_normal_lpdf(datum, mu, sigma));
    }
    return result;
}


template<class Hypers> 
void HierarchyNNW<Hypers>::draw(){
    Eigen::MatrixXd tau_new = stan::math::wishart_rng( this->hypers->get_nu(),
        this->hypers->get_tau0(), this->rng );
    Eigen::MatrixXd sigma = this->state[1].inverse();
    EigenRowVec mu_new = stan::math::multi_normal_rng( this->hypers->get_mu0(),
        sigma*(1/this->hypers->get_lambda()), this->rng );

     this->state[0] = mu_new;
     this->state[1] = tau_new;
}


template<class Hypers> 
Eigen::VectorXd HierarchyNNW<Hypers>::eval_marg(const Eigen::MatrixXd &data){
    unsigned int n = data.rows();
    Eigen::VectorXd result(n);
    unsigned int dim = data.cols();
    double nu = this->hypers->get_nu();
    double lambda = this->hypers->get_lambda();

    EigenRowVec mu_n = this->hypers->get_mu0();
    double nu_n = 2*nu - dim + 1;
    Eigen::MatrixXd sigma_n = this->hypers->get_tau0().inverse() *
        ( nu-(dim-1)/2 ) * lambda/(lambda+1);
    EigenRowVec row;
    for(int i = 0; i < n; i++){
        // use multi_student_t_lpdf(datum, nu, mu, Sigma)
        row = data.row(i);
        result(i) = exp( stan::math::multi_student_t_lpdf(row, nu_n,
            mu_n, sigma_n) );
    }
    return result;
}


template<class Hypers> 
std::vector<Eigen::MatrixXd> HierarchyNNW<Hypers>::normal_wishart_update(
    const Eigen::MatrixXd &data, const EigenRowVec &mu0, const double lambda,
    const Eigen::MatrixXd &tau0, const double nu){
    unsigned int n = data.rows();
    Eigen::MatrixXd lambda_post(1,1), nu_post(1,1);
    EigenRowVec mubar = data.colwise().mean();
    lambda_post(0,0) = lambda + n;
    nu_post(0,0) = nu + 0.5 * n;
    EigenRowVec mu_post = (lambda*mu0 + n*mubar) * (1/(lambda+n));

    // Compute tau_post
    Eigen::MatrixXd tau_temp = Eigen::MatrixXd::Zero(n, n);
    for(unsigned int i = 0; i < n; i++){
        EigenRowVec datum = data.row(i);
        tau_temp += (datum-mubar).transpose()*(datum-mubar); // column-times-row
    }
    tau_temp += (n*lambda/(n+lambda)) * (mubar-mu0).transpose()*(mubar-mu0);
    tau_temp = 0.5*tau_temp + tau0.inverse();
    Eigen::MatrixXd tau_post = tau_temp.inverse();

    return std::vector<Eigen::MatrixXd>{mu_post,lambda_post,tau_post,nu_post};
}


template<class Hypers> 
void HierarchyNNW<Hypers>::sample_given_data(const Eigen::MatrixXd &data){
    // Get current values of parameters
    EigenRowVec mu0 = this->hypers->get_mu0();
    double lambda = this->hypers->get_lambda();
    Eigen::MatrixXd tau0 = this->hypers->get_tau0();
    double nu = this->hypers->get_nu();
    std::vector<Eigen::MatrixXd> temp = normal_wishart_update(data, mu0, lambda,
        tau0, nu);
    EigenRowVec mu_post = temp[0];
    double lambda_post = temp[1](0,0);
    Eigen::MatrixXd tau_post = temp[2];
    double nu_post = temp[3](0,0);
    // Get a sample
    Eigen::MatrixXd tau_new = stan::math::wishart_rng(nu_post, tau_post,
        this->rng);
    Eigen::MatrixXd tau_inv = tau_new.inverse();
    EigenRowVec mu_new = stan::math::multi_normal_rng(mu_post,
        tau_inv*(1/lambda_post), this->rng);
    
    this->state[0] = mu_new;
    this->state[1] = tau_new;
}


#endif // HIERARCHYNNW_IMP_HPP
