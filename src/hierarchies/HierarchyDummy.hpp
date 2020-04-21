#ifndef HIERARCHYDUMMY_HPP
#define HIERARCHYDUMMY_HPP

#include "HierarchyBase.hpp"

// Dummmy multivariate hierarchy for testing purposes
//
// state  = mu, sigma (vector, matrix)
// hypers = mu0, lambda0 (vector, matrix)

template<class Hypers>
class HierarchyDummy : public HierarchyBase<Hypers> {
protected:
    std::vector<Eigen::MatrixXd> dummy_update(const Eigen::MatrixXd &data,
        const Eigen::VectorXd &mu0, const Eigen::MatrixXd &lambda0);

public:
    bool is_multivariate() const override {return true;}

    // Destructor and constructor
    ~HierarchyDummy() = default;

    HierarchyDummy(std::shared_ptr<Hypers> hypers_) {
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

#include "HierarchyDummy.imp.hpp"

#endif // HIERARCHYDUMMY_HPP
