#ifndef HIERARCHYNW_IMP_HPP
#define HIERARCHYNW_IMP_HPP

#include "HierarchyNW.hpp"

// TODO check everything

template<class Hypers> 
void HierarchyNW<Hypers>::set_tau_and_utilities(const Eigen::MatrixXd &tau){
    if(this->state.size() == 1){ // e.g. if the hierarchy is being initialized
        this->state.push_back(tau);
    }
    else {
        this->state[1] = tau;
    }

    chol_factor = LLT<Eigen::MatrixXd>(tau);
    chol_factor_eval = chol_factor.matrixL();
    Eigen::VectorXd diag = chol_factor_eval.diagonal();
    log_det = 2 * log(diag.array()).sum();
}


template<class Hypers> 
Eigen::VectorXd HierarchyNW<Hypers>::like(const Eigen::MatrixXd &data){
    // instead of using the inefficient stan::math::multi_normal_lpdf()
    int n = data.cols();
    double out = tau_log_det * n;

    Eigen::MatrixXd chol_sigma = tau_chol_factor_eval;

    std::vector<double> loglikes(n);
    for(int i = 0; i < n; i++){
        loglikes[i] = (chol_sigma * (x[i] - mu)).squaredNorm();
    }

    out -= std::accumulate(loglikes.begin(), loglikes.end(), 0.0);

    //return 0.5 * out; // TODO ?????
    return Eigen::MatrixXd(n);

    //Eigen::VectorXd result(data.rows());
    //for(int i = 0; i < data.rows(); i++){
    //    result(i) = exp( stan::math::multi_normal_lpdf(data(i),this->state[0],
    //        inverse) );
    //}
    //return result;
}


template<class Hypers> 
void HierarchyNW<Hypers>::draw(){
    this->state[0] = stan::math::multi_normal_rng(
        this->hypers->get_mu0(),
        (this->state->get_tau0_inv()) * (1/this->hypers->get_lambda()),
        this->rng
    );

    Eigen::MatrixXd sigma_new = stan::math::wishart_rng( this->hypers->get_nu(),
        this->hypers->get_tau0_inv(), this->rng );
    set_tau_and_utilities(sigma_new.inverse());
}


template<class Hypers> 
Eigen::VectorXd HierarchyNW<Hypers>::eval_marg(const Eigen::MatrixXd &data){
    Eigen::VectorXd result(data.cols());
    // TODO to do lol
    return result;
}


template<class Hypers> 
std::vector<Eigen::MatrixXd> HierarchyNW<Hypers>::normal_wishart_update(
    const Eigen::MatrixXd &data, const Eigen::VectorXd &mu0,
    const double lambda, const Eigen::MatrixXd &tau0, const double nu){
    int n = data.cols();
    Eigen::MatrixXd lambda_post(1,1), nu_post(1,1);
    Eigen::VectorXd mubar = data.rowwise().mean();

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
    set_tau_and_utilities(sigma_new.inverse());
}


#endif // HIERARCHYNW_IMP_HPP
