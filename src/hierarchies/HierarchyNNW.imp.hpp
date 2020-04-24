#ifndef HIERARCHYNNW_IMP_HPP
#define HIERARCHYNNW_IMP_HPP

#include "HierarchyNNW.hpp"

// TODO check everything

template<class Hypers> 
void HierarchyNNW<Hypers>::set_tau_and_utilities(const Eigen::MatrixXd &tau){
    if(this->state.size() == 1){ // e.g. if the hierarchy is being initialized
        this->state.push_back(tau);
    }
    else {
        this->state[1] = tau;
    }

    tau_chol_factor = Eigen::LLT<Eigen::MatrixXd>(tau);
    tau_chol_factor_eval = tau_chol_factor.matrixL();
    Eigen::VectorXd diag = tau_chol_factor_eval.diagonal();
    tau_log_det = 2 * log(diag.array()).sum();
}


template<class Hypers> 
Eigen::VectorXd HierarchyNNW<Hypers>::like(const Eigen::MatrixXd &data){
    // instead of using the inefficient stan::math::multi_normal_lpdf()
    int n = data.rows();
    Eigen::VectorXd result(n);
    for(unsigned int i = 0; i < n; i++){
        result(i) = 0.5 * std::exp( n*tau_log_det * (tau_chol_factor_eval*
            (data.row(i) - this->state[0])).squaredNorm() );
    }
    return result;
}


template<class Hypers> 
void HierarchyNNW<Hypers>::draw(){
    Eigen::MatrixXd tau_new = stan::math::wishart_rng( this->hypers->get_nu(),
        this->hypers->get_tau0(), this->rng );
    Eigen::MatrixXd mu_new = stan::math::multi_normal_rng(
        this->hypers->get_mu0(),
        this->state[1].inverse()*(1/this->hypers->get_lambda()),
        this->rng );

     this->state[0] = mu_new;
     set_tau_and_utilities(tau_new);
}


template<class Hypers> 
Eigen::VectorXd HierarchyNNW<Hypers>::eval_marg(const Eigen::MatrixXd &data){
    unsigned int n = data.rows();
    Eigen::VectorXd result(n);
    unsigned int dim = data.cols();
    double nu = this->hypers->get_nu();
    double lambda = this->hypers->get_lambda();

    Eigen::VectorXd mu_n = this->hypers->get_mu0();
    double nu_n = 2*nu - dim + 1;
    Eigen::MatrixXd sigma_n = this->hypers->get_tau0().inverse() *
        ( nu-(dim-1)/2 ) * lambda/(lambda+1);

    for(int i = 0; i < n; i++){
        // use multi_student_t_lpdf(datum, nu, mu, Sigma)
        result(i) = exp( stan::math::multi_student_t_lpdf(data.row(i), nu_n,
            mu_n, sigma_n) );
    }
    return result;
}


template<class Hypers> 
std::vector<Eigen::MatrixXd> HierarchyNNW<Hypers>::normal_wishart_update(
    const Eigen::MatrixXd &data, const Eigen::VectorXd &mu0,
    const double lambda, const Eigen::MatrixXd &tau0, const double nu){
    unsigned int n = data.rows();
    Eigen::MatrixXd lambda_post(1,1), nu_post(1,1);
    Eigen::VectorXd mubar = data.colwise().mean();

    lambda_post(0,0) = lambda + n;
    nu_post(0,0) = nu + n;
    Eigen::VectorXd mu_post = (lambda*mu0 + n*mubar) * (1/lambda+n);

    // Compute tau_post
    Eigen::MatrixXd tau_temp = tau0.inverse();
    for(unsigned int i = 0; i < n; i++){
        tau_temp += (data.row(i)-mubar)*(data.row(i)-mubar).transpose();
    }
    tau_temp += (nu*lambda/(nu+lambda)) * (mubar-mu0)*(mubar-mu0).transpose();
    Eigen::MatrixXd tau_post = tau_temp.inverse();
    
    return std::vector<Eigen::MatrixXd>{mu_post,lambda_post,tau_post,nu_post};
}


template<class Hypers> 
void HierarchyNNW<Hypers>::sample_given_data(const Eigen::MatrixXd &data){
    // Get current values of parameters
    Eigen::VectorXd mu0 = this->hypers->get_mu0();
    double lambda = this->hypers->get_lambda();
    Eigen::MatrixXd tau0 = this->hypers->get_tau0();
    double nu = this->hypers->get_nu();

    std::vector<Eigen::MatrixXd> temp = normal_wishart_update(data, mu0, lambda,
        tau0, nu);

    Eigen::VectorXd mu_post = temp[0];
    double lambda_post = temp[1](0,0);
    Eigen::MatrixXd tau_post = temp[2];
    double nu_post = temp[3](0,0);

    // Get a sample
    Eigen::MatrixXd tau_new = stan::math::wishart_rng(nu_post, tau_post,
        this->rng);
    Eigen::VectorXd mu_new = stan::math::multi_normal_rng(mu_post,
        tau_new.inverse()*(1/lambda_post), this->rng);

    this->state[0] = mu_new;
    set_tau_and_utilities(tau_new);
}


#endif // HIERARCHYNNW_IMP_HPP
