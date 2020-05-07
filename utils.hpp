#ifndef UTILS_HPP
#define UTILS_HPP

#include <Eigen/Dense>
#include <fstream>

void fill_eigen_matrix_from_file(Eigen::Ref<Eigen::MatrixXd> mat,
    const std::string &filename){
    // Needs space-separated values!
    std::ifstream istr(filename);

    if(istr.is_open())
    {
        for (int i = 0; i < mat.rows(); i++)
            for (int j = 0; j < mat.cols(); j++)
            {
                double val;
                istr >> val;
                mat(i,j) = val;
            }
        istr.close();
    }
    // TODO mat2 move into mat?
}

#endif // UTILS_HPP