#ifndef UTILS_HPP
#define UTILS_HPP

#include <fstream>
#include <numeric>
#include <vector>
#include <random>
#include <Eigen/Dense>
#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/util/delimited_message_util.h>
#include <stan/math/prim/mat.hpp>
#include "precs/precmat.hpp"
#include "../protos/cpp/state.pb.h"

using namespace Eigen;

void delete_row(MatrixXd* x, int ind);

void delete_elem(VectorXd* x, int ind);

MatrixXd delete_row(const MatrixXd &x, int ind);

VectorXd delete_elem(const VectorXd &x, int ind);

MatrixXd vstack(const std::vector<VectorXd>& rows);

template <typename T>
T loadTextProto(std::string filename)
{
    std::ifstream ifs(filename);
    google::protobuf::io::IstreamInputStream iis(&ifs);
    T out;
    auto success = google::protobuf::TextFormat::Parse(&iis, &out);
    if (!success)
        std::cout << "An error occurred in 'loadTextProto'; success: " 
        << success << std::endl;
    return out;
}

double o_multi_normal_prec_lpdf(
    const VectorXd &x, const VectorXd &mu, const PrecMat &sigma);

double o_multi_normal_prec_lpdf(
    const std::vector<VectorXd> &x, const VectorXd &mu, const PrecMat &sigma);

double trunc_normal_rng(
    double mu, double sigma, double lower, double upper,
    std::mt19937_64 &rng);

double trunc_normal_lpdf(
    double x, double mu, double sigma, double lower, double upper);

void to_proto(const MatrixXd &mat, EigenMatrix *out);

void to_proto(const VectorXd &vec, EigenVector *out);

std::vector<VectorXd> to_vector_of_vectors(const MatrixXd& mat);

MatrixXd pairwise_dist_sq(const MatrixXd &x, const MatrixXd &y);

MatrixXd pairwise_dist_sq(const MatrixXd &x);

#endif 