#include "simulate_straus.hpp"

Eigen::MatrixXd simulate_strauss_moller(
    const Eigen::MatrixXd &ranges, double beta, double gamma, double R) {

    gsl_rng *r;
    r = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(r, rand() % 100000);

    // Initialise point pattern
    Point2Pattern ExamplePattern(
        (ranges(0, 0) - 1) * 2,
        (ranges(1, 0) + 1) * 2,
        (ranges(0, 1) - 1) * 2,
        (ranges(1, 1) + 1) * 2, 9, 9);

    Sampler ExampleSimulator;

    StraussProcess ExampleProcess(
        ranges(0, 0), ranges(1, 0), // range of x
        ranges(0, 1), ranges(1, 1), // range of y
        beta, gamma, R);

    // Generate perfect sample of Strauss process
    ExampleSimulator.Sim(&ExamplePattern, &ExampleProcess, r);

    return ExamplePattern.to_eigen_mat();
}