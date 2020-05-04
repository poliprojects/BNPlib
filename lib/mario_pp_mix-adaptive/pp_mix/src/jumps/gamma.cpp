#include "gamma.hpp"

double GammaJump::sample_tilted(double u)
{
    // double s, v;
    // while (1 > 0) {
    //     s = gamma_rng(alpha, beta, Rng::Instance().get());
    //     v = uniform_rng(0, 1, Rng::Instance().get());
    //     if (std::log(v) < -s * u)
    //         return s;
    // }
    return gamma_rng(alpha, beta + u, Rng::Instance().get());
}

double GammaJump::sample_given_data(
    int ndata, double curr, double u)
{

    double out;
    int nh = ndata;
    double prop = curr + uniform_rng(-0.1, 0.1, Rng::Instance().get());

    double num = std::log(prop) * nh - u * prop +
                 gamma_lpdf(prop, alpha, beta);

    double den = std::log(curr) * nh - u * curr +
                 gamma_lpdf(curr, alpha, beta);

    if (std::log(uniform_rng(0, 1, Rng::Instance().get())) < num - den) 
        out = prop;
    else
        out = curr;
    return out;
}

double GammaJump::laplace(double u) {
    return std::pow(beta, alpha) / std::pow(beta + u, alpha);
}