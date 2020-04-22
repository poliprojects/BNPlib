#ifndef HIERARCHYNW_HPP
#define HIERARCHYNW_HPP

#include "HierarchyBase.hpp"

// Normal-Wishart multivariate distribution
//
// state  = mu, tau (vector, matrix)
// hypers = mu0, lambda, tau0, nu (vector, scalar, matrix, scalar)

template<class Hypers>
class HierarchyNW : public HierarchyBase<Hypers> {
protected:
    std::vector<Eigen::MatrixXd> nw_update(const Eigen::MatrixXd &data,
        const Eigen::VectorXd &???);
    Eigen::MatrixXd inverse; // TODO

public:
    bool is_multivariate() const override {return true;}

    // Destructor and constructor
    ~HierarchyNW() = default;

    HierarchyNW(std::shared_ptr<Hypers> hypers_) {
        this->hypers = hypers_;
        dim = this->hypers.get_mu0().size();
        this->state.push_back( this->hypers.get_mu0() );
        this->state.push_back( this->hypers.get_lambda() *
        	Eigen::Matrix<double, dim, dim>::Identity() );
    }

    Eigen::VectorXd eval_marg(const Eigen::MatrixXd &data) override;
    Eigen::VectorXd like(const Eigen::MatrixXd &data) override;

    void draw() override;

    void sample_given_data(const Eigen::MatrixXd &data) override;
};

#include "HierarchyNW.imp.hpp"

#endif // HIERARCHYNW_HPP
