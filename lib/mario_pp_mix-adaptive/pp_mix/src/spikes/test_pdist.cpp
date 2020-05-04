#include <Eigen/Dense>
#include <iostream>
using namespace Eigen;


MatrixXd pairwise_dist_sq(const MatrixXd &x, const MatrixXd &y)
{
    MatrixXd D(x.rows(), y.rows());
    int i=0;
    for (int i = 0; i < y.rows(); i++)
        D.col(i) = (x.rowwise() - y.row(i)).rowwise().squaredNorm();

    return D;
}

int main() {
    MatrixXd x = Eigen::MatrixXd::Random(5, 3);
    std::cout << "x: \n" << x << std::endl;
    MatrixXd y = Eigen::MatrixXd::Random(10, 3);
    std::cout << "y: \n" << y << std::endl;

    MatrixXd pdist1 = pairwise_dist_sq(x, y);

    std::cout << "pdist1\n" << pdist1 << std::endl;

    MatrixXd pdist2(5, 10);
    for (int i=0; i < 5; i++) {
        for (int j=0; j < 10; j++) {
            pdist2(i, j) = (x.row(i) - y.row(j)).array().square().sum();
        }
    }

    std::cout << "pdist2\n" << pdist2 << std::endl;

    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 10; j++)
        {   
            std::cout << "i: " << i << ", j: " << j; 
            assert(std::abs(pdist1(i, j) - pdist2(i, j)) < 1e-5);
            std::cout << " pass" << std::endl;
        }
    }
}