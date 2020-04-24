#ifndef INCLUDES_MAIN_HPP
#define INCLUDES_MAIN_HPP

#include "src/algorithms/Neal2.hpp"
#include "src/algorithms/Neal8.hpp"
#include "src/hierarchies/HierarchyNNIG.hpp"
#include "src/hierarchies/HierarchyNNW.hpp"
#include "src/hyperparameters/HypersFixedNNIG.hpp"
#include "src/hyperparameters/HypersFixedNNW.hpp"
#include "src/mixtures/DirichletMixture.hpp"
#include "src/mixtures/PitYorMixture.hpp"
#include "src/runtime/Factory.hpp"


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
}


#endif // INCLUDES_MAIN_HPP
