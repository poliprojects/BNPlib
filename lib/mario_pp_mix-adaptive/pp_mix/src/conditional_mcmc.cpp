#include "conditional_mcmc.hpp"
#include "conditional_mcmc_imp.hpp"

MultivariateConditionalMCMC::MultivariateConditionalMCMC(
    BasePP *pp_mix, BaseJump *h, BasePrec *g) : 
        ConditionalMCMC<BaseMultiPrec, PrecMat, VectorXd>()
{
    set_pp_mix(pp_mix);
    set_jump(h);
    set_prec(dynamic_cast<BaseMultiPrec*>(g));
}

void MultivariateConditionalMCMC::get_state_as_proto(
        google::protobuf::Message *out_)
{
    using namespace google::protobuf::internal;

    MultivariateMixtureState *out = down_cast<MultivariateMixtureState *>(out_);
    out->set_ma(a_means.rows());
    out->set_mna(na_means.rows());
    out->set_mtot(a_means.rows() + na_means.rows());

    for (int i = 0; i < a_means.rows(); i++)
    {
        EigenVector *mean;
        EigenMatrix *prec;
        mean = out->add_a_means();
        prec = out->add_a_precs();

        to_proto(a_means.row(i).transpose(), mean);
        to_proto(a_precs[i].get_prec(), prec);
    }

    for (int i = 0; i < na_means.rows(); i++)
    {
        EigenVector *mean;
        EigenMatrix *prec;
        mean = out->add_na_means();
        prec = out->add_na_precs();

        to_proto(na_means.row(i).transpose(), mean);
        to_proto(na_precs[i].get_prec(), prec);
    }

    to_proto(a_jumps, out->mutable_a_jumps());
    to_proto(na_jumps, out->mutable_na_jumps());

    *out->mutable_clus_alloc() = {clus_alloc.data(), clus_alloc.data() + ndata};

    out->set_u(u);

    PPState pp_params;
    pp_mix->get_state_as_proto(&pp_params);
    out->mutable_pp_state()->CopyFrom(pp_params);
}

UnivariateConditionalMCMC::UnivariateConditionalMCMC(
    BasePP *pp_mix, BaseJump *h, BasePrec *g) : 
        ConditionalMCMC<BaseUnivPrec, double, double>()
{
    set_pp_mix(pp_mix);
    set_jump(h);
    set_prec(dynamic_cast<BaseUnivPrec *>(g));
}

void UnivariateConditionalMCMC::get_state_as_proto(
            google::protobuf::Message *out_) 
{
    using namespace google::protobuf::internal;

    UnivariateMixtureState *out = down_cast<UnivariateMixtureState *>(out_);
    out->set_ma(a_means.rows());
    out->set_mna(na_means.rows());
    out->set_mtot(a_means.rows() + na_means.rows());

    to_proto(Map<VectorXd>(a_means.data(), a_means.rows()), 
             out->mutable_a_means());
    to_proto(Map<VectorXd>(na_means.data(), na_means.rows()), 
             out->mutable_na_means());

    EigenVector* precs = out->mutable_a_precs();
    precs->set_size(a_precs.size());
    *precs->mutable_data() = {a_precs.begin(), a_precs.end()};

    precs = out->mutable_na_precs();
    precs->set_size(na_precs.size());
    *precs->mutable_data() = {na_precs.begin(), na_precs.end()};

    to_proto(a_jumps, out->mutable_a_jumps());
    to_proto(na_jumps, out->mutable_na_jumps());

    *out->mutable_clus_alloc() = {clus_alloc.data(), clus_alloc.data() + ndata};

    out->set_u(u);

    PPState pp_params;
    pp_mix->get_state_as_proto(&pp_params);
    out->mutable_pp_state()->CopyFrom(pp_params);
}