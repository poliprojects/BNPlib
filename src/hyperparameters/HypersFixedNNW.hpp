#ifndef HYPERSFIXEDNNW_HPP
#define HYPERSFIXEDNNW_HPP

#include <cassert>
#include <Eigen/Dense>

class HypersFixedNNW {
private:
    using EigenRowVec = Eigen::Matrix<double, 1, Eigen::Dynamic>;

    EigenRowVec mu0;
    double lambda;
    Eigen::MatrixXd tau0;
    double nu;

    void check_state_validity(){
        unsigned int dim = mu0.size();
        assert(lambda > 0);
        assert(dim == tau0.rows());
        assert(tau0.rows() == tau0.cols());
        assert(nu > dim-1);

        // Check if tau0 is symmetric positive semi definite
        assert( tau0.isApprox(tau0.transpose()) );
        Eigen::LLT<Eigen::MatrixXd> llt(tau0);
        assert( llt.info() != Eigen::NumericalIssue );
    }

public:
    // Destructor and constructor
    ~HypersFixedNNW() = default;
    HypersFixedNNW() = default;
    HypersFixedNNW(const EigenRowVec &mu0_, const double lambda_,
        const Eigen::MatrixXd &tau0_, const double nu_):
        mu0(mu0_), lambda(lambda_), tau0(tau0_), nu(nu_) {
        
        check_state_validity();
        }

    // Getters
    const EigenRowVec get_mu0(){return mu0;}
    const double get_lambda(){return lambda;}
    const Eigen::MatrixXd get_tau0(){return tau0;}
    const double get_nu(){return nu;}

    // Setters
    void set_mu0(const EigenRowVec &mu0_){
    	assert(mu0_.size() == mu0.size());
    	mu0 = mu0_;
    }
    void set_lambda(const double lambda_){
    	assert(lambda_ > 0);
    	lambda = lambda_;
    }
    void set_tau0(const Eigen::MatrixXd &tau0_) {
    	assert(tau0_.rows() == tau0_.cols());
    	assert(mu0.size() == tau0_.rows());
        assert( tau0.isApprox(tau0.transpose()) );
        Eigen::LLT<Eigen::MatrixXd> llt(tau0);
        assert( llt.info() != Eigen::NumericalIssue );
    	tau0 = tau0_;
    }
    void set_nu(const double nu_){
    	assert(nu_ > mu0.size()-1);
    	nu = nu_;
    }
};


#endif // HYPERSFIXEDNNW_HPP
