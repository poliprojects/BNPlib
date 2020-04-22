#ifndef HIERARCHYNW_HPP
#define HIERARCHYNW_HPP

#include "HierarchyBase.hpp"

// Normal-Wishart multivariate distribution
//
// state  = mu, Tau (vector, matrix)
// hypers = mu0, lambda, Tau0, nu (vector, scalar, matrix, scalar)

template<class Hypers>
class HierarchyNW : public HierarchyBase<Hypers> {
protected:
    std::vector<Eigen::MatrixXd> nw_update(const Eigen::MatrixXd &data,
        const Eigen::VectorXd &mu0, const Eigen::MatrixXd &lambda0);

public:
    bool is_multivariate() const override {return true;}

    // Destructor and constructor
    ~HierarchyNW() = default;

    HierarchyNW(std::shared_ptr<Hypers> hypers_) {
        this->hypers = hypers_;
        Eigen::VectorXd mu(3);
        mu << 2.9, 2.9, 2.9;
        this->state.push_back(mu);
        Eigen::MatrixXd sig(3,3);
        sig << 1.0, 0.0, 0.0,
               0.0, 1.0, 0.0,
               0.0, 0.0, 1.0;
        this->state.push_back(sig);
    }

    Eigen::VectorXd eval_marg(const Eigen::MatrixXd &datum) override;
    Eigen::VectorXd like(const Eigen::MatrixXd &datum) override;

    void draw() override;

    void sample_given_data(const Eigen::MatrixXd &data) override;
};

#include "HierarchyNW.imp.hpp"

#endif // HIERARCHYNW_HPP
