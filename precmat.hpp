#ifndef COVMAT_HPP
#define COVMAT_HPP

#include <iostream>
#include <Eigen/Dense>
using namespace Eigen;


/*
 * This class stores a precision matrix, along with some functionals
 * of it used to compute the pdf of a normal r.v.
 */
class PrecMat
{
protected:
    MatrixXd prec;
    LLT<MatrixXd> cho_factor;
    MatrixXd cho_factor_eval;
    double log_det;

public:
    PrecMat() {}
    ~PrecMat() {}

    PrecMat(const MatrixXd &prec);

    MatrixXd get_prec() const;

    LLT<MatrixXd> get_cho_factor() const;

    const MatrixXd &get_cho_factor_eval() const;

    double get_log_det() const;
};

std::ostream &operator<<(std::ostream &output, const PrecMat &p);

#endif
