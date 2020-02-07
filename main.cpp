#include <iostream>

#include <boost/random/random_number_generator.hpp>
#include <boost/random/detail/qrng_base.hpp>

#include "includes_main.hpp"

int main() {
    double mean1=5;
    double mean2=5;
    double sd1=1;
    double sd2=1;
    std::mt19937 rng_base;
    std::vector<double> data(40);
    int half = data.size()/2;

    for(int i = 0; i < half; i++){
        data[i]      = stan::math::normal_rng(mean1, sd1, rng_base);
        data[i+half] = stan::math::normal_rng(mean2, sd2, rng_base);
    }

    HypersFixed hy(4,1.5,2,2);
    SimpleMixture mix(1.0);
    //Neal2<NNIGHierarchy, HypersFixed, SimpleMixture> sampler2(data, mix, hy);
    Neal8<NNIGHierarchy, HypersFixed, SimpleMixture> sampler8(data, 3, mix, hy);
    //std::cout << "Running Neal2" << std::endl;
    //sampler2.run();
    std::cout << "Running Neal8" << std::endl;
    sampler8.run();
    //sampler8.write_clustering_to_file();
    std::vector<double> grid = {1,1.5,2,2.5,3,3.5,4};
    //sampler8.eval_density(grid);
    //sampler8.write_density_to_file();

    //for(int i = 0; i < 40; i++){
      //  std::cout << i << ": " << data[i] << std::endl;
    //}

    return 0;

}
