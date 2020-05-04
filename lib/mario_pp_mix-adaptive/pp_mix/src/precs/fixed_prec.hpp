#ifndef FIXED_COV_HPP
#define FIXED_COV_HPP

#include "base_prec.hpp"

class FixedUnivPrec : public BaseUnivPrec
{
protected:
   double sigma;

public:
   FixedUnivPrec(double sigma) : sigma(sigma) {}

   ~FixedUnivPrec() {}

   double sample_prior() override { return sigma; }

   double sample_given_data(
       const std::vector<double> &data, const double &curr,
       const VectorXd &mean) override { return sigma; }
};

class FixedPrec : public BaseMultiPrec
{
protected:
   int dim;
   double sigma;

public:
   FixedPrec(int dim, double sigma) : dim(dim), sigma(sigma) {}

   ~FixedPrec() {}

   PrecMat sample_prior() override;

   PrecMat sample_given_data(
       const std::vector<VectorXd> &data, const PrecMat &curr,
       const VectorXd &mean) override;
};


#endif 