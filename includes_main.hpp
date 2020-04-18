#ifndef INCLUDES_MAIN_HPP
#define INCLUDES_MAIN_HPP

#include "src/hyperparameters/HypersFixedNNIG.hpp"
#include "src/hyperparameters/HypersDummy.hpp"
#include "src/algorithms/Neal2.hpp"
#include "src/algorithms/Neal8.hpp"
#include "src/algorithms/Factory.hpp"
#include "src/hierarchies/HierarchyNNIG.hpp"
#include "src/hierarchies/HierarchyDummy.hpp"
#include "src/mixtures/DirichletMixture.hpp"

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
