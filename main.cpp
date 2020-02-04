#include <iostream>
#include "Neal8_NNIG.hpp"
#include "Neal2_NNIG.hpp"
#include "SimpleMixture.hpp"
#include "HypersFixed.hpp"
#include "NNIGHierarchy.hpp"
#include <vector>
#include <Eigen/Dense> 
#include <stan/math/prim/mat.hpp>
#include <boost/random/random_number_generator.hpp>
#include <boost/random/detail/qrng_base.hpp>
#include "includes_main.hpp"

int main() {
    double mean1=5.5;
    double mean2=2.5;
    double sd1=1;
    double sd2=1;
    std::mt19937 rng_base;
    std::vector<double> data(40);
    int half = data.size()/2;

    for (int i=0; i<half; i++) {
        data[i]      = stan::math::normal_rng(mean1, sd1, rng_base);
        data[i+half] = stan::math::normal_rng(mean2, sd2, rng_base);
    }

    HypersFixed hy(4,1.5,2,2);
    SimpleMixture mix(1.0);

    Neal8<NNIGHierarchy, HypersFixed, SimpleMixture> sampler(data, 3, mix, hy);
    sampler.run();

    for(int i=0; i<40; i++)
        std::cout << i << ": " << data[i] << std::endl;

    return 0;

}

