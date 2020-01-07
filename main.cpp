#include <iostream>
#include "Neal8_NNIG.hpp"
#include <vector>
#include <armadillo>
#include <Eigen/Dense> 
#include <stan/math/prim/mat.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/random/detail/qrng_base.hpp>

#include "includes.hpp"

int main() {

  double means=5;
  double sds=1;
  std::mt19937 rng_base;
  std::vector<double> data(20);
  for (int i=0; i< data.size(); i++) {

  data[i] = stan::math::normal_rng(means, sds, rng_base);
  //std::cout<<data[i]<<std::endl;
}


    HypersFixed hy(4,1,1,1);
    SimpleMixture mix(5.0);

    Neal8<NNIGHierarchy, HypersFixed, SimpleMixture> sampler(data, 3, mix, hy);
    sampler.run();

	std::cout << "Test" << std::endl;
	return 0;

}
