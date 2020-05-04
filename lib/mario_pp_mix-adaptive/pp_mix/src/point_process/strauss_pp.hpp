#ifndef STRAUSS_PP
#define STRAUSS_PP

#include "base_pp.hpp"
#include "../simulate_straus.hpp"
#include "../adaptive_metropolis.hpp"
#include "../../protos/cpp/params.pb.h"
#include <google/protobuf/stubs/casts.h>
#include <boost/math/distributions/chi_squared.hpp>

class StraussPP: public BasePP {
 protected:
    double beta, gamma, R;
    bool fixed_params = false;
    StraussParams::Priors priors;
    AdaptiveMetropolis<double, double> am_beta;

    double sqrt_chisq_quantile;


 public:

    StraussPP() {}

    StraussPP(StraussParams::Priors priors);

    StraussPP(double beta, double gamma, double R);

    StraussPP(
        StraussParams::Priors priors, double beta, double gamma, double R);

    ~StraussPP() {}

    void initialize() override;

    double dens(const MatrixXd &x, bool log = true) override;

    double dens_from_pdist(const MatrixXd& dists, double beta_, double gamma_,
                           double R_, bool log=true);

    double papangelou(
        MatrixXd xi, const MatrixXd &x, bool log = true) override;

    VectorXd phi_star_rng() override;

    double phi_star_dens(VectorXd xi, bool log = true) override;

    void update_hypers(const MatrixXd &active, const MatrixXd &non_active) override;

    void get_state_as_proto(google::protobuf::Message *out);

    double estimate_mean_proposal_sigma();
};

#endif