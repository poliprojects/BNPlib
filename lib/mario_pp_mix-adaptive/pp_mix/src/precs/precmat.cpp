#include "precmat.hpp"


PrecMat::PrecMat(const MatrixXd &prec): prec(prec) {
    cho_factor = LLT<MatrixXd>(prec);
    cho_factor_eval = cho_factor.matrixL();
    const VectorXd& diag = cho_factor_eval.diagonal();
    log_det = 2 * log(diag.array()).sum();
}

MatrixXd PrecMat::get_prec() const
{
    return prec;
}

LLT<MatrixXd> PrecMat::get_cho_factor() const
{
    return cho_factor;
}

const MatrixXd& PrecMat::get_cho_factor_eval() const
{
    return cho_factor_eval;
}

double PrecMat::get_log_det() const
{
    return log_det;
}