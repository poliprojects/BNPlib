#ifndef HIERARCHYNNIG_HPP
#define HIERARCHYNNIG_HPP

#include "HierarchyBase.hpp"


//! Normal Normal-InverseGamma hierarchy for univariate data.

//! This class represents a hiearchy, i.e. a cluster, whose data are distributed
//! according to a normal likelihood, the parameters of which have a Normal-In-
//! verseGamma centering distribution. That is:
//! f(xi|mu,sig^2) = N(mu,sig^2)
//!     (mu,sig^2) ~ G
//!              G ~ MM(M, G0)
//!             G0 = N-IG
//! Therefore the state of this hierarchy is (mu, sigma) and their hyperparame-
//! ters contained in the Hypers object are (mu_0, lambda, alpha, beta).
//! Note that this hierarchy is conjugate, thus the marginal and the posterior
//! distribution are available in closed form.

//! \param Hypers Name of the hyperparameters class

template<class Hypers>
class HierarchyNNIG : public HierarchyBase<Hypers> {
protected:
    using HierarchyBase<Hypers>::state;
    using HierarchyBase<Hypers>::hypers;

    // AUXILIARY TOOLS
    //! Raises error if the state values are not valid w.r.t. their own domain
    void check_state_validity() override;

    //! ...
    std::vector<double> normal_gamma_update(const Eigen::VectorXd &data,
        const double mu0, const double alpha0, const double beta0,
        const double lambda);

public:
    //! Returns true if the hierarchy models multivariate data (here, false)
    bool is_multivariate() const override {return false;}

    // DESTRUCTOR AND CONSTRUCTORS
    ~HierarchyNNIG() = default;
    HierarchyNNIG() = default;
    HierarchyNNIG(std::shared_ptr<Hypers> hypers_) {
        hypers = hypers_;
        state = std::vector<Eigen::MatrixXd>(2,Eigen::MatrixXd(1,1));
        state[0](0,0) = hypers->get_mu0();
        state[1](0,0) = 1;
    }

    // EVALUATION FUNCTIONS
    //! Evaluates the likelihood of data in the given points
    Eigen::VectorXd like(const Eigen::MatrixXd &data) override;
    //! Evaluates the marginal distribution of data in the given points
    Eigen::VectorXd eval_marg(const Eigen::MatrixXd &data) override;

    // SAMPLING FUNCTIONS
    //! Generates new values for state from its prior distribution
    void draw() override;
    //! Generates new values for state from its posterior distribution
    void sample_given_data(const Eigen::MatrixXd &data) override;
};


#include "HierarchyNNIG.imp.hpp"

#endif // HIERARCHYNNIG_HPP
