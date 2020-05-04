#include "../covs/covmat.hpp"
#include "../utils.hpp"

#include <vector>
#include <stan/math/prim/mat.hpp>

int main() {
    MatrixXd Sigma_(2,2);
    Sigma_ << 1, 0.3, 0.3, 1; 

    CovMat sigma(Sigma_);

    VectorXd x(2);
    x << 1, 5;

    VectorXd mu(2);
    mu << 0.5, 3;

    std::cout << "ours: " << multi_normal_lpdf(x, mu, sigma) << std::endl;
    std::cout << "stan: " << stan::math::multi_normal_lpdf(x, mu, Sigma_) << std::endl;

    double out1 = 0;
    double out2;

    std::vector<VectorXd> x_vec;

    for (int i=0; i < 20; i++) {
        out1 += multi_normal_lpdf(x, mu, sigma);
        x_vec.push_back(x);
    }

    out2 = multi_normal_lpdf(x_vec, mu, sigma);

    std::cout << "out1: " << out1 << std::endl;
    std::cout << "out2: " << out2 << std::endl;

}