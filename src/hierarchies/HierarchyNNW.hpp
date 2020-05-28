#ifndef HIERARCHYNNW_HPP
#define HIERARCHYNNW_HPP

#include "HierarchyBase.hpp"


//! ...
// Normal-Wishart multivariate distribution
//
// state  = mu, tau (vector, matrix)
// hypers = mu0, lambda, tau0, nu (vector, scalar, matrix, scalar)

//! \param Hypers Name of the hyperparameters class

template<class Hypers>
class HierarchyNNW : public HierarchyBase<Hypers> {
protected:
    using HierarchyBase<Hypers>::state;
    using HierarchyBase<Hypers>::hypers;
    using EigenRowVec = Eigen::Matrix<double, 1, Eigen::Dynamic>;

    // UTILITIES FOR LIKELIHOOD COMPUTATION
    //! ...
    Eigen::LLT<Eigen::MatrixXd> tau_chol_factor;
    //! ...
    Eigen::MatrixXd tau_chol_factor_eval;
    //! ...
    double tau_log_det;

    // AUXILIARY TOOLS
    //! Raises error if parameters are not valid w.r.t. their own domain
    void check_state_validity() override;
    //! Special setter for the second part of the 
    void set_tau_and_utilities(const Eigen::MatrixXd &tau);

    //! ...
    std::vector<Eigen::MatrixXd> normal_wishart_update(
        const Eigen::MatrixXd &data, const EigenRowVec &mu0,
        const double lambda, const Eigen::MatrixXd &tau0, const double nu);

public:
    //! Returns true if the hierarchy models multivariate data (here, true)
    bool is_multivariate() const override {return true;}

    // DESTRUCTOR AND CONSTRUCTORS
    ~HierarchyNNW() = default;
    HierarchyNNW() = default;
    HierarchyNNW(std::shared_ptr<Hypers> hypers_) {
        hypers = hypers_;
        unsigned int dim = hypers->get_mu0().size();
        state.push_back( hypers->get_mu0() );
        set_tau_and_utilities( hypers->get_lambda() *
            Eigen::MatrixXd::Identity(dim, dim) );
    }

    // EVALUATION FUNCTIONS
    //! Evaluates model likelihood in the given points
    Eigen::VectorXd like(const Eigen::MatrixXd &data) override;
    //! Evaluates marginal distribution of data in the given points
    Eigen::VectorXd eval_marg(const Eigen::MatrixXd &data) override;

    // SAMPLING FUNCTIONS
    //! ...
    void draw() override;
    //! ...
    void sample_given_data(const Eigen::MatrixXd &data) override;

    // GETTERS AND SETTERS
    //! \param state_ State value to set
    //! \param check  If true, a state validity check occurs after assignment
    void set_state(const std::vector<Eigen::MatrixXd> &state_,
    	bool check = true) override {
        state[0] = state_[0];
        set_tau_and_utilities(state_[1]);
        if(check){
        	check_state_validity();
        }
    }
};

#include "HierarchyNNW.imp.hpp"

#endif // HIERARCHYNNW_HPP
