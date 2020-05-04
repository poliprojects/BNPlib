#ifndef NREP_PP_HPP
#define NREP_PP_HPP

#include <Eigen/Dense>
#include <boost/math/distributions/gamma.hpp>
#include "base_pp.hpp"


class NrepPP: public BasePP {
protected:
    double tau, u, p;

public:
    NrepPP() {}
    ~NrepPP() {}

    NrepPP(double u, double p);

    void initialize() override;

    /*
    * This method performs the calibration of \tau proposed in
    * Section 3.1.1 in Quinlan et al. (2017)
    */
    void calibrate();

    double dens(const MatrixXd &x, bool log = true);

    double papangelou(
        MatrixXd xi, const MatrixXd &x, bool log = true);

    VectorXd phi_star_rng();

    double phi_star_dens(VectorXd xi, bool log = true);

    void update_hypers(const MatrixXd &active, const MatrixXd &non_active);

    void get_state_as_proto(google::protobuf::Message *out);

    /*
    * this method is called by ConditionalMCMC::sample_means()
    * it should return a value for the standard deviation of the proposal
    * density in such a way that it gives sufficient high probability
    * to regions where the repulsion is lower.
    */
    double estimate_mean_proposal_sigma();


    double multi_trunc_normal_lpdf(const VectorXd& x);
};

#endif