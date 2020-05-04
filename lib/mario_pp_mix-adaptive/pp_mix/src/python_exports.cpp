#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <Eigen/Dense>
#include <deque>
#include <string>

#include "conditional_mcmc.hpp"
#include "factory.hpp"
#include "utils.hpp"
#include "../protos/cpp/params.pb.h"


namespace py = pybind11;

std::deque<py::bytes> run_pp_mix_univ(
    int burnin, int niter, int thin, const Eigen::MatrixXd &data,
    Params params, const Eigen::MatrixXd &ranges)
{

    int log_every = 200;
    std::deque<py::bytes> out;
    BasePP *pp_mix = make_pp(params);
    BaseJump *h = make_jump(params);
    BasePrec *g = make_prec(params);
    pp_mix->set_ranges(ranges);

    UnivariateConditionalMCMC sampler(pp_mix, h, g);
    const std::vector<double> datavec(data.data(), data.data() + data.size());
    sampler.initialize(datavec);

    for (int i = 0; i < burnin; i++)
    {
        sampler.run_one();
        if ((i + 1) % log_every == 0)
        {
            py::print("Burnin, iter #", i + 1, " / ", burnin);
        }
    }

    for (int i = 0; i < niter; i++)
    {
        sampler.run_one();
        if (i % thin == 0)
        {
            std::string s;
            UnivariateMixtureState curr;
            sampler.get_state_as_proto(&curr);
            curr.SerializeToString(&s);
            out.push_back((py::bytes)s);
        }

        if ((i + 1) % log_every == 0)
        {
            py::print("Running, iter #", i + 1, " / ", niter);
        }
    }
    return out;
}

std::deque<py::bytes> run_pp_mix_multi(
    int burnin, int niter, int thin, const Eigen::MatrixXd &data,
    Params params, const Eigen::MatrixXd &ranges)
{

    int log_every = 200;
    std::deque<py::bytes> out;
    BasePP *pp_mix = make_pp(params);
    BaseJump *h = make_jump(params);
    BasePrec *g = make_prec(params);
    pp_mix->set_ranges(ranges);

    MultivariateConditionalMCMC sampler(pp_mix, h, g);
    std::vector<Eigen::VectorXd> datavec = to_vector_of_vectors(data);
    sampler.initialize(datavec);

    for (int i = 0; i < burnin; i++)
    {
        sampler.run_one();
        if ((i + 1) % log_every == 0)
        {
            py::print("Burnin, iter #", i + 1, " / ", burnin);
        }
    }

    for (int i = 0; i < niter; i++)
    {
        sampler.run_one();
        if (i % thin == 0)
        {
            std::string s;
            MultivariateMixtureState curr;
            sampler.get_state_as_proto(&curr);
            curr.SerializeToString(&s);
            out.push_back((py::bytes)s);
        }

        if ((i + 1) % log_every == 0)
        {
            py::print("Running, iter #", i + 1, " / ", niter);
        }
    }
    return out;
}

std::deque<py::bytes> _run_pp_mix(
        int burnin, int niter, int thin, const Eigen::MatrixXd& data,
        std::string serialized_params) {

    Eigen::MatrixXd ranges(2, data.cols());
    ranges.row(0) = data.colwise().minCoeff();
    ranges.row(1) = data.colwise().maxCoeff();
    ranges *= 2;
    std::cout << "ranges: \n" << ranges << std::endl;

    Params params;
    params.ParseFromString(serialized_params);

    if (data.rows() == 1 || data.cols() == 1)
        return run_pp_mix_univ(burnin, niter, thin, data, params, ranges);
    else
        return run_pp_mix_multi(burnin, niter, thin, data, params, ranges);
}


Eigen::MatrixXd _simulate_strauss2D(
    const Eigen::MatrixXd &ranges, double beta, double gamma, double R)
{
    return simulate_strauss_moller(ranges, beta, gamma, R);
}

PYBIND11_MODULE(pp_mix_cpp, m)
{
    m.doc() = "aaa"; // optional module docstring

    m.def("_run_pp_mix", &_run_pp_mix, "aaa");

    m.def("_simulate_strauss2D", &_simulate_strauss2D, "aaa");
}
