#include <pybind11/pybind11.h>

#include "../../includes.hpp"
#include "run_NNIG_Dir.cpp"
#include "run_NNW_Dir.cpp"
#include "estimates_NNIG_Dir.cpp"
#include "estimates_NNIG_Dir_grid.cpp"



PYBIND11_MODULE(bnplib, m)
{
    m.doc() = "TODO docstring";
    m.def("run_NNIG_Dir", &run_NNIG_Dir, "TODO docstring");
    m.def("run_NNW_Dir", &run_NNW_Dir, "TODO docstring");
    m.def("estimates_NNIG_Dir", &estimates_NNIG_Dir, "TODO docstring");
    m.def("estimates_NNIG_Dir_grid", &estimates_NNIG_Dir_grid, "TODO docstr");
}
