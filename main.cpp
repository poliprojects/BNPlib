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

  double mean1=5;
  double mean2=0;
  double sd1=1;
  double sd2=1;
  std::mt19937 rng_base;
  std::vector<double> data(40);
  int half = data.size()/2;

  for (int i=0; i< half; i++) {

  data[i]      = stan::math::normal_rng(mean1, sd1, rng_base);
  data[i+half] = stan::math::normal_rng(mean2, sd2, rng_base);
  
  //std::cout<<data[i]<<std::endl;
}
	//for(auto &c : data)
	//	std::cout << c << " ";
	//std::cout << std::endl;
	//return 0;

    HypersFixed hy(4,1,1,1);
    SimpleMixture mix(1.0);

    Neal8<NNIGHierarchy, HypersFixed, SimpleMixture> sampler(data, 3, mix, hy);
    sampler.run();

	//std::cout << "Test" << std::endl;
	return 0;

}
