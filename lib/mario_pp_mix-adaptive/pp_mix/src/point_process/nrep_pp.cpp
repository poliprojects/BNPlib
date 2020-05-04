#include "nrep_pp.hpp"


NrepPP::NrepPP(double u, double p): u(u), p(p) { }

void NrepPP::initialize()
{
    c_star = 1.0;
    calibrate();
}

void NrepPP::calibrate() 
{
    boost::math::gamma_distribution<double> gamma((1.0 * dim) / 2, 0.5);
    double q = quantile(complement(gamma, 1.0 - p));
    tau = q / (-log(1.0 - u));
    std::cout << "tau: " << tau << std::endl;
}

double NrepPP::dens(const MatrixXd &x, bool log)
{   

    double out = 0.0;
    if ((x.size() == 1 && dim == 1) || (x.rows() == 1 & dim > 1) ||
        (x.cols() == 1 && dim > 1))
    {   
        out = multi_trunc_normal_lpdf(VectorXd::Map(x.data(), x.size()));
    }
    else
    {
        for (int i=0; i < x.rows(); i++)
            out += multi_trunc_normal_lpdf(x.row(i).transpose());

        MatrixXd pdist = pairwise_dist_sq(x);

        out += (ArrayXXd::Constant(pdist.rows(), pdist.rows(), 1.0) -
                (-0.5 / tau * pdist.array()).exp())
                   .log()
                   .sum();
    }

    if (! log)
        out = exp(out);

    return out;
}

double NrepPP::papangelou(MatrixXd xi, const MatrixXd &x, bool log) 
{
    double out = 0.0;
    if (xi.cols() == 1)
        xi.transposeInPlace();

    for (int i=0; i < xi.rows(); i++)
        out += multi_trunc_normal_lpdf(xi.row(i).transpose());

    MatrixXd dists = pairwise_dist_sq(xi, x);
    ArrayXXd id = ArrayXXd::Constant(dists.rows(), dists.cols(), 1.0);
    out += (id - (-0.5 / tau * dists.array()).exp()).log().sum();
    if (!log)
        out = exp(out);
    
    return out;
}

VectorXd NrepPP::phi_star_rng()
{
    VectorXd out(dim);
    for (int i=0; i < dim; i++) {
        out(i) = trunc_normal_rng(0.0, 1.0, ranges(0, i), ranges(1, i),
                                  Rng::Instance().get());
    }

    return out;
}

double NrepPP::phi_star_dens(VectorXd xi, bool log) 
{

    double out = multi_trunc_normal_lpdf(xi);
    if (!log)
        out = exp(out);
    return out;
}

void NrepPP::update_hypers(const MatrixXd &active, const MatrixXd &non_active)
{
    return;
}

void NrepPP::get_state_as_proto(google::protobuf::Message *out) 
{
    return;
}

double NrepPP::estimate_mean_proposal_sigma() {
    return 0.25;
}

double NrepPP::multi_trunc_normal_lpdf(const VectorXd &x) {
    double out = 0.0;

    for (int i=0; i < dim; i++)
        out += trunc_normal_lpdf(x(i), 0.0, 1.0, ranges(0, i), ranges(1, i));

    return out;
}