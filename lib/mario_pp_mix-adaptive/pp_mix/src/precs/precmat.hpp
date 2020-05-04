#ifndef COVMAT_HPP
#define COVMAT_HPP

#include <Eigen/Dense>
using namespace Eigen;


class PrecMat {
 protected:
   MatrixXd prec;
   LLT<MatrixXd> cho_factor;
   MatrixXd cho_factor_eval;
   double log_det;
   double univariate_val;
   bool is_univariate = false;

 public:
   PrecMat() {}
   ~PrecMat() {}

   PrecMat(const MatrixXd &prec);

   PrecMat(const double &prec): univariate_val(prec) {
      is_univariate = true;
   }

   MatrixXd get_prec() const;

   LLT<MatrixXd> get_cho_factor() const;

   const MatrixXd& get_cho_factor_eval() const;

   double get_log_det() const;
};

#endif 