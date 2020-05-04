#include "factory.hpp"

BasePP *make_pp(const Params &params)
{
    BasePP *out;
    if (params.has_strauss())
        out = make_strauss(params.strauss());
    else if (params.has_nrep())
        out = make_nrep(params.nrep());

    return out;
}

BasePP *make_strauss(const StraussParams &params) {
    BasePP* out;

    if (params.fixed_params())
        out = new StraussPP(
            params.init().beta(), params.init().gamma(), params.init().r());
    else if (params.has_init() and params.has_prior())
        out = new StraussPP(
            params.prior(), params.init().beta(), params.init().gamma(),
            params.init().r());
    else
        out = new StraussPP(params.prior());

    return out;
}

BasePP *make_nrep(const NrepParams &params)
{
    return new NrepPP(params.u(), params.p());
}

BaseJump *make_jump(const Params &params)
{
    BaseJump *out;
    if (params.has_gamma_jump())
    {
        out = make_gamma_jump(params.gamma_jump());
    }
    return out;
}

BaseJump *make_gamma_jump(const GammaParams &params) {
    return new GammaJump(params.alpha(), params.beta());
}

BasePrec *make_prec(const Params &params)
{
    BasePrec *out;
    if (params.has_fixed_multi_prec())
        out = make_fixed_prec(params.fixed_multi_prec());
    else if (params.has_wishart())
        out = make_wishart(params.wishart());
    else if (params.has_fixed_univ_prec())
        out = make_fixed_prec(params.fixed_univ_prec());
    else if (params.has_gamma_prec())
        out = make_gamma_prec(params.gamma_prec());

    return out;
}

BasePrec *make_fixed_prec(const FixedMultiPrecParams &params)
{
    return new FixedPrec(params.dim(), params.sigma());
}

BasePrec *make_wishart(const WishartParams &params) {
    double sigma = 1.0;
    if (!params.sigma__case())
        sigma = params.sigma();

    return new Wishart(params.nu(), params.dim(), sigma);
}

BasePrec *make_fixed_prec(const FixedUnivPrecParams &params)
{
    return new FixedUnivPrec(params.sigma());
}

BasePrec *make_gamma_prec(const GammaParams &params) {
    return new GammaPrec(params.alpha(), params.beta());
}
