#include "utils.hpp"

void delete_row(MatrixXd *x, int ind)
{
    int nrow = x->rows() - 1;
    int dim = x->cols();
    if (ind < nrow)
    {
        x->block(ind, 0, nrow - ind, dim) =
            x->block(ind + 1, 0, nrow - ind, dim);
    }

    x->conservativeResize(nrow, dim);
}

void delete_elem(VectorXd *x, int ind)
{
    int size = x->size() - 1;
    if (ind < size)
        x->segment(ind, size - ind) = x->segment(ind + 1, size - ind);

    x->conservativeResize(size);
}


MatrixXd delete_row(const MatrixXd &x, int ind) {
    MatrixXd out = x;
    int nrow = x.rows() - 1;
    int dim = x.cols();
    if (ind < nrow)
    {
        out.block(ind, 0, nrow - ind, dim) =
            out.block(ind + 1, 0, nrow - ind, dim);
    }

    out.conservativeResize(nrow, dim);
    return out;
}

VectorXd delete_elem(const VectorXd &x, int ind)
{   
    VectorXd out = x;
    int size = x.size() - 1;
    if (ind < size)
        out.segment(ind, size - ind) = out.segment(ind + 1, size - ind);

    out.conservativeResize(size);
    return out;
}

MatrixXd vstack(const std::vector<VectorXd> &rows) {
    int nrows = rows.size();
    int ncols = rows[0].size();

    MatrixXd out(nrows, ncols);
    for (int i=0; i < nrows; i++)
        out.row(i) = rows[i].transpose();

    return out;
}


double o_multi_normal_prec_lpdf(
    const VectorXd &x, const VectorXd &mu, const PrecMat &sigma)
{
    // std::cout << "sigma.cho_factor:\n " << sigma.get_cho_factor_eval() << std::endl;
    double out = sigma.get_log_det();
    out -= (sigma.get_cho_factor_eval() * (x - mu)).squaredNorm();
    return 0.5 * out;
}

double o_multi_normal_prec_lpdf(
    const std::vector<VectorXd> &x, const VectorXd &mu, const PrecMat &sigma)
{
    int n = x.size();
    double out = sigma.get_log_det() * n;

    const MatrixXd& cho_sigma = sigma.get_cho_factor_eval();

    std::vector<double> loglikes(n);
    for (int i = 0; i < n; i++)
    {
        loglikes[i] = (cho_sigma * (x[i] - mu)).squaredNorm();
    }

    out -= std::accumulate(loglikes.begin(), loglikes.end(), 0.0);

    return 0.5 * out;
}

// generate from truncated normal by rejection sampling
// !! might not be the best idea
double trunc_normal_rng(
    double mu, double sigma, double lower, double upper,
    std::mt19937_64 &rng)
{
    while (true)
    {
        double val = stan::math::normal_rng(mu, sigma, rng);
        if (val <= upper && val >= lower)
            return val;
    }
}

double trunc_normal_lpdf(double x, double mu, double sigma, double lower, double upper)
{
    if ((x < lower) || (x > upper))
        return stan::math::NEGATIVE_INFTY;

    double out = stan::math::normal_lpdf(x, mu, sigma);
    out -= stan::math::log_diff_exp(
        stan::math::normal_lcdf(upper, mu, sigma),
        stan::math::normal_lcdf(lower, mu, sigma));

    return out;
}

void to_proto(const MatrixXd &mat, EigenMatrix* out) {
    out->set_rows(mat.rows());
    out->set_cols(mat.cols());
    *out->mutable_data() = {mat.data(), mat.data() + mat.size()};
}

void to_proto(const VectorXd &vec, EigenVector* out)
{
    out->set_size(vec.size());
    *out->mutable_data() = {vec.data(), vec.data() + vec.size()};
}

std::vector<VectorXd> to_vector_of_vectors(const MatrixXd &mat)
{
    std::vector<VectorXd> out(mat.rows());
    for (int i = 0; i < mat.rows(); i++)
        out[i] = mat.row(i).transpose();

    return out;
}

MatrixXd pairwise_dist_sq(const MatrixXd &x, const MatrixXd &y)
{
    MatrixXd D(x.rows(), y.rows());
    int i = 0;
    for (int i = 0; i < y.rows(); i++)
        D.col(i) = (x.rowwise() - y.row(i)).rowwise().squaredNorm();

    return D;
}

MatrixXd pairwise_dist_sq(const MatrixXd &x)
{
    return pairwise_dist_sq(x, x);
}