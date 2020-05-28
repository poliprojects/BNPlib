#ifndef HIERARCHYNNIG_HPP
#define HIERARCHYNNIG_HPP

#include "HierarchyBase.hpp"


//! ...
// Normal likelihoood, Normal-InverseGamma hierarchy, that is:
// f ~ N(mu,sig^2)
// (mu,sig^2) ~ G
// G ~ DP(M, G0)  with G0 = N-IG
//
// state  = mu, sigma
// hypers = mu_0, lambda, alpha, beta

//! \param Hypers Name of the hyperparameters class

template<class Hypers>
class HierarchyNNIG : public HierarchyBase<Hypers> {
protected:
    using HierarchyBase<Hypers>::state;
    using HierarchyBase<Hypers>::hypers;

    // AUXILIARY TOOLS
    //! Raises error if parameters are not valid w.r.t. their own domain
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
    //! Evaluates model likelihood in the given points
    Eigen::VectorXd like(const Eigen::MatrixXd &data) override;
    //! Evaluates marginal distribution of data in the given points
    Eigen::VectorXd eval_marg(const Eigen::MatrixXd &data) override;

    // SAMPLING FUNCTIONS
    //! ...
    void draw() override;
    //! ...
    void sample_given_data(const Eigen::MatrixXd &data) override;
};


#include "HierarchyNNIG.imp.hpp"

#endif // HIERARCHYNNIG_HPP
