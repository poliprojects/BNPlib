#include "strauss_pp.hpp"

StraussPP::StraussPP(StraussParams::Priors priors) : priors(priors)
{

    beta = (priors.beta_u() + priors.beta_l()) / 2.0;
    beta = (priors.gamma_u() + priors.gamma_l()) / 2.0;
    beta = (priors.r_u() + priors.r_l()) / 2.0;
    fixed_params = false;
}

StraussPP::StraussPP(
        double beta, double gamma, double R):
            beta(beta), gamma(gamma), R(R) {
    fixed_params = true;
}

StraussPP::StraussPP(
    StraussParams::Priors priors, double beta, double gamma, double R) : 
        priors(priors), beta(beta), gamma(gamma), R(R)
{
    fixed_params = false;
}

void StraussPP::initialize() {
    am_beta = AdaptiveMetropolis<double, double>(1);
    am_beta.init();
    c_star = beta;
    boost::math::chi_squared chisq(dim);
    sqrt_chisq_quantile = sqrt(quantile(complement(chisq, 0.1)));
}

double StraussPP::dens(const MatrixXd &x, bool log) {
    double out;
    // check if it's jut one point
    if ((x.size() == 1 && dim == 1) || (x.rows() == 1 & dim > 1) ||
        (x.cols() == 1 && dim > 1))
        out = std::log(beta);
    else {
        int npoints = x.rows();
        out = std::log(beta) * npoints;
        MatrixXd pdist = pairwise_dist_sq(x);

        out += std::log(gamma) * (pdist.array() < R*R).count();
    }
    if (! log)
        out = std::exp(out);
    return out;
}

double StraussPP::dens_from_pdist(const MatrixXd &dists, double beta_, double gamma_,
                                  double R_, bool log)
{
    double out;
    if (dists.rows() == 1)
        out = beta_;
    else  {
        int npoints = dists.rows();
        double out = std::log(beta_) * npoints;
        out += std::log(gamma_) * (dists.array() < R_ * R_).count();
    }
    if (!log)
        out = std::exp(out);
    return out;
}

double StraussPP::papangelou(
    MatrixXd xi, const MatrixXd &x, bool log)
{   
    double out;
    if (xi.cols() == 1)
        xi.transposeInPlace();
    MatrixXd dists = pairwise_dist_sq(xi, x);
    out = std::log(gamma) * (dists.array() < R * R).count();
    if (!log)
        out = std::exp(out);
    return out;
}

VectorXd StraussPP::phi_star_rng() {
    VectorXd out(dim);
    for (int i=0; i < dim; i++) {
        out(i) = uniform_rng(ranges(0, i), ranges(1, i), Rng::Instance().get());
    }
    return out;
}

double StraussPP::phi_star_dens(VectorXd xi, bool log)
{
    double out = beta;
    if (log)
        out = std::log(out);

    return out;
}


double StraussPP::estimate_mean_proposal_sigma() {
    return R / sqrt_chisq_quantile;
}

void StraussPP::update_hypers(
        const MatrixXd& active, const MatrixXd& non_active) {
    if (fixed_params)
        return;

    if (active.cols() != 2)
        return;

    // std::cout << "StraussPP::update_hypers" << std::endl;
    int Ma = active.rows();
    int Mna = non_active.rows();
    int dim = active.cols();
    MatrixXd all_points(Ma + Mna, dim);
    all_points.block(0, 0, Ma, dim) = active;
    all_points.block(Ma, 0, Mna, dim) = non_active;

    MatrixXd dists = pairwise_dist_sq(all_points);
    MatrixXd aux_var, aux_dists;

    double scale;
    double prop, arate, lower, upper;
    
    // UPDATE BETA
    upper = priors.beta_u();
    lower = priors.beta_l();
    scale = (upper - lower) / 5;

    prop = trunc_normal_rng(beta, scale, lower, upper, Rng::Instance().get());

    double proprate = trunc_normal_lpdf(prop, beta, scale, lower, upper) -
                      trunc_normal_lpdf(beta, prop, scale, lower, upper);

    arate = proprate;

    double likrate = dens_from_pdist(dists, prop, gamma, R) -
                     dens_from_pdist(dists, beta, gamma, R);
    arate += likrate;

    aux_var = simulate_strauss_moller(ranges, prop, gamma, R);
    aux_dists = pairwise_dist_sq(aux_var);

    double aux_lik_rate = dens_from_pdist(aux_dists, beta, gamma, R) -
                          dens_from_pdist(aux_dists, prop, gamma, R);
    arate += aux_lik_rate;
    
    if (log(uniform_rng(0, 1, Rng::Instance().get())) < arate)
    {
        beta = prop;
    }
    
    // std::cout << "####################################\n"
    //           << std::endl;
    // // UPDATE Gamma
    // upper = priors.gamma_u();
    // lower = priors.gamma_l();
    // scale = (upper - lower) / 3;
    // prop = trunc_normal_rng(gamma, scale, lower, upper, Rng::Instance().get());

    // arate = trunc_normal_lpdf(prop, gamma, scale, lower, upper) -
    //         trunc_normal_lpdf(gamma, prop, scale, lower, upper);

    // arate += dens_from_pdist(dists, beta, prop, R) -
    //          dens_from_pdist(dists, beta, gamma, R);

    // aux_var = simulate_strauss_moller(ranges, beta, prop, R);
    // aux_dists = pairwise_dist_sq(aux_var);

    // arate += dens_from_pdist(aux_dists, beta, gamma, R) -
    //          dens_from_pdist(aux_dists, beta, prop, R);

    // if (uniform_rng(0, 1, Rng::Instance().get()) < arate) {
    //     gamma = prop;
    //     std::cout << "accepted gamma" << std::endl;
    // }

    // // UPDATE R
    // upper = priors.r_u();
    // lower = priors.r_l();
    // scale = (upper - lower) / 3;
    // prop = trunc_normal_rng(R, scale, lower, upper, Rng::Instance().get());

    // arate = trunc_normal_lpdf(prop, R, scale, lower, upper) -
    //         trunc_normal_lpdf(R, prop, scale, lower, upper);

    // arate += dens_from_pdist(dists, beta, gamma, prop) -
    //          dens_from_pdist(dists, beta, gamma, R);

    // aux_var = simulate_strauss_moller(ranges, beta, gamma, prop);
    // aux_dists = pairwise_dist_sq(aux_var);

    // arate += dens_from_pdist(aux_dists, beta, gamma, R) -
    //          dens_from_pdist(aux_dists, beta, gamma, prop);

    // if (uniform_rng(0, 1, Rng::Instance().get()) < arate) {
    //     R = prop;
    //     std::cout << "accepted R" << std::endl;
    // }
}

void StraussPP::get_state_as_proto(google::protobuf::Message *out)
{
    StraussState state;
    state.set_beta(beta);
    state.set_gamma(gamma);
    state.set_r(R);
    state.set_birth_prob(birth_prob);
    state.set_birth_arate(birth_arate);

    using namespace google::protobuf::internal;
    down_cast<PPState *>(out)->mutable_strauss_state()->CopyFrom(state);
}