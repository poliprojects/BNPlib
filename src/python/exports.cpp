#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include "../../includes.hpp"
#include "run_NNIG_Dir.cpp"
#include "run_NNW_Dir.cpp"
#include "estimates_NNIG_Dir.cpp"


PYBIND11_MODULE(bnplib, m)
{
    m.doc() = "TODO docstring";
    m.def("run_NNIG_Dir", &run_NNIG_Dir, "TODO docstring");
    m.def("run_NNW_Dir", &run_NNW_Dir, "TODO docstring");
    m.def("estimates_NNIG_Dir", &estimates_NNIG_Dir, "TODO docstring");
}
