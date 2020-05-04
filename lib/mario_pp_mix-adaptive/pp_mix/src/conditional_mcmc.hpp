#ifndef CONDITIONAL_MCMC
#define CONDITIONAL_MCMC

// #include <omp.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <functional>

#include <Eigen/Dense>
#include <stan/math/prim/mat.hpp>
#include <google/protobuf/message.h>

#include "rng.hpp"
#include "point_process/base_pp.hpp"
#include "jumps/base_jump.hpp"
#include "precs/base_prec.hpp"
#include "precs/precmat.hpp"
#include "utils.hpp"
#include "../protos/cpp/state.pb.h"

using namespace Eigen;


template<class Prec, typename prec_t, typename data_t>
class ConditionalMCMC {
 protected:
     int dim;
     int ndata;
     std::vector<data_t> data;
     double prior_ratio, lik_ratio;
     std::vector<std::vector<data_t>> data_by_clus;

     // STATE
     int nclus;
     double u;
     VectorXi clus_alloc;
     VectorXd a_jumps, na_jumps;
     MatrixXd a_means, na_means;
     std::vector<prec_t> a_precs, na_precs;

     // DISTRIBUTIONS
     BasePP *pp_mix;
     BaseJump *h;
     Prec *g;

     // FOR DEBUGGING
     bool verbose = false;
     int acc_mean = 0;
     int tot_mean = 0;

 public:
     ConditionalMCMC() {}
     ~ConditionalMCMC()
     {
         delete pp_mix;
         delete h;
         delete g;
    }

    ConditionalMCMC(BasePP * pp_mix, BaseJump * h, Prec * g);

    void set_pp_mix(BasePP* pp_mix) {this->pp_mix = pp_mix;}
    void set_jump(BaseJump* h) {this->h = h;}
    void set_prec(Prec* g) {this->g = g;}

    void initialize(const std::vector<data_t> &data);

    void run_one();

    void sample_allocations_and_relabel();

    void sample_means();

    void sample_vars();

    void sample_jumps();

    virtual void get_state_as_proto(google::protobuf::Message *out_) = 0;

    void print_debug_string();

    void _relabel();

    void set_verbose() { verbose = !verbose; }

    virtual double normal_lpdf_single(
        const data_t &x, const VectorXd &mu, const prec_t &sigma)  = 0;

    virtual double normal_lpdf_single_multi(
        const std::vector<data_t> &x, const VectorXd &mu, 
        const prec_t &sigma) = 0;

    virtual void set_dim(const data_t& datum) = 0;

    double mean_acceptance_rate() {
        return (1.0 * acc_mean) / (1.0 * tot_mean);
    }
};


class MultivariateConditionalMCMC: public ConditionalMCMC<
    BaseMultiPrec, PrecMat, VectorXd> {

public:
    MultivariateConditionalMCMC() {}

    MultivariateConditionalMCMC(BasePP *pp_mix, BaseJump *h, BasePrec *g);

    void get_state_as_proto(google::protobuf::Message *out_) override;

    double normal_lpdf_single(
        const VectorXd &x, const VectorXd &mu, const PrecMat &sigma)
    {
        return o_multi_normal_prec_lpdf(x, mu, sigma);
    }

    double normal_lpdf_single_multi(
        const std::vector<VectorXd> &x, const VectorXd &mu, const PrecMat &sigma)
    {
        return o_multi_normal_prec_lpdf(x, mu, sigma);
    }

    void set_dim(const VectorXd& datum) {
        dim = datum.size();
    }

};

class UnivariateConditionalMCMC : public ConditionalMCMC<
    BaseUnivPrec, double, double>
{
public:
    UnivariateConditionalMCMC() {}

    UnivariateConditionalMCMC(BasePP *pp_mix, BaseJump *h, BasePrec *g);

    void get_state_as_proto(google::protobuf::Message *out_) override;

    double normal_lpdf_single(
        const double &x, const VectorXd &mu, const double &sigma)
    {
        return stan::math::normal_lpdf(x, mu(0), 1.0 / sigma);
    }

    double normal_lpdf_single_multi(
        const std::vector<double> &x, const VectorXd &mu, const double &sigma)
    {
        return stan::math::normal_lpdf(x, mu(0), 1.0 / sigma);
    }

    void set_dim(const double &datum)
    {
        dim = 1;
    }
};

#include "conditional_mcmc_imp.hpp"

#endif