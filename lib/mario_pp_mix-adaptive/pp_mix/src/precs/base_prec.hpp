#ifndef BASE_COV_HPP
#define BASE_COV_HPP

#include <vector>
#include <Eigen/Dense>

#include "precmat.hpp"

using namespace Eigen;

class BasePrec
{
public:
    virtual ~BasePrec(){};
};

class BaseUnivPrec: public BasePrec {
public:
    virtual ~BaseUnivPrec() {};
    
    virtual double sample_prior() = 0;

    virtual double sample_given_data(
        const std::vector<double> &data, const double &curr,
        const VectorXd &mean) = 0;
};


class BaseMultiPrec: public BasePrec {
public:
    virtual ~BaseMultiPrec(){};

    virtual PrecMat sample_prior() = 0;

    virtual PrecMat sample_given_data(
        const std::vector<VectorXd> &data, const PrecMat &curr,
        const VectorXd &mean) = 0;
};

#endif