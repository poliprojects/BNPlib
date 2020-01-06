#include <iostream>
#include "Neal8_NNIG.hpp"
#include <vector>
#include <armadillo>
#include <Eigen/Dense> 
#include <stan/math/prim/mat.hpp>

#include "includes.hpp"

int main() {

  arma::vec means(1);
  arma::vec sds(1);
  means(1)=5;
  sds(1)=1;
  
  //boost::random_number_generator rng_t;

  std::vector<double> data(20);
  //for (int i=0; i< data.size(); i++) {
  //data[i] = stan::math::normal_rng(means(1), sds(1), rng_t);}


    //HypersFixed hy();
    //SimpleMixture mix(5.0);
    //Neal8<NNIGHierarchy<HypersFixed>, SimpleMixture> sampler(hy,mix);

	std::cout << "Test" << std::endl;
	return 0;

}
