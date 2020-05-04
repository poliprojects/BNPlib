#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
//#include <pybind11/eigen.h>

#include "../../includes.hpp"
#include "run_NNIG.cpp"


PYBIND11_MODULE(bnplib, m)
{
    m.doc() = "TODO docstring";
    m.def("run_NNIG", &run_NNIG, "TODO docstring");
}
