#include "gamma.hpp"


GammaPrec::GammaPrec(double alpha, double beta): alpha(alpha), beta(beta) {}


double GammaPrec::sample_prior() {
    return stan::math::gamma_rng(alpha, beta, Rng::Instance().get());
}

double GammaPrec::sample_given_data(
    const std::vector<double> &data, const double &curr,
    const VectorXd &mean) {
    
    double m = mean(0);
    std::vector<double> x_min_mu;
    std::transform(data.begin(), data.end(), std::back_inserter(x_min_mu),
                   [m](const double x) {return (x-m) * (x-m);});

    double sum_squares = std::accumulate(x_min_mu.begin(), x_min_mu.end(), 0.0);

    return stan::math::gamma_rng(
        alpha + (1.0 * data.size()) / 2,
        beta + 0.5 * sum_squares,
        Rng::Instance().get());
    return 1;
}