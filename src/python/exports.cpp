#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include "../../includes.hpp"
#include "run_NNIG_Dir.cpp"
#include "run_NNW_Dir.cpp"
#include "estimates_NNIG_Dir.cpp"
#include "estimates_NNW_Dir.cpp"

//! \file

//! Builds the Python version of the library.

namespace py = pybind11;

PYBIND11_MODULE(bnplibpy, m)
{
    m.doc() = "Nonparametric library for cluster and density estimation";

    m.def("run_NNIG_Dir", &run_NNIG_Dir,
        "Runs algorithm for an NNIG + Dirichlet mixture model",
        py::arg("mu0"), py::arg("lambda_"), py::arg("alpha0"), py::arg("beta0"),
        py::arg("totalmass"), py::arg("datafile"), py::arg("algo"),
        py::arg("collfile"), py::arg("init"), py::arg("rng"), py::arg("maxit"),
        py::arg("burn"), py::arg("n_aux"));

    m.def("run_NNW_Dir", &run_NNW_Dir,
        "Runs algorithm for an NNW + Dirichlet mixture model",
        py::arg("mu0"), py::arg("lambda_"), py::arg("tau0"), py::arg("nu"),
        py::arg("totalmass"), py::arg("datafile"), py::arg("algo"),
        py::arg("collfile"), py::arg("init"), py::arg("rng"), py::arg("maxit"),
        py::arg("burn"), py::arg("n_aux"));

    m.def("estimates_NNIG_Dir", &estimates_NNIG_Dir,
        "Cluster and density estimates for an NNIG + Dirichlet mixture model",
        py::arg("mu0"), py::arg("lambda_"), py::arg("alpha0"), py::arg("beta0"),
        py::arg("totalmass"), py::arg("grid"), py::arg("algo"),
        py::arg("collfile"), py::arg("densfile"), py::arg("clustfile"),
        py::arg("only"));

    m.def("estimates_NNW_Dir", &estimates_NNW_Dir,
        "Cluster and density estimates for an NNW + Dirichlet mixture model",
        py::arg("mu0"), py::arg("lambda_"), py::arg("tau0"), py::arg("nu"),
        py::arg("totalmass"), py::arg("grid"), py::arg("algo"),
        py::arg("collfile"), py::arg("densfile"), py::arg("clustfile"),
        py::arg("only"));
}
