#include <iostream>

#include <boost/random/random_number_generator.hpp>
#include <boost/random/detail/qrng_base.hpp>

#include "includes_main.hpp"

int main(){
    double mean1 = 5.0;
    double mean2 = 5.0;
    double sd1   = 1.0;
    double sd2   = 1.0;
    std::mt19937 rng_base;
    std::vector<double> data(40);
    int half = data.size()/2;

    for(int i = 0; i < half; i++){
        data[i]      = stan::math::normal_rng(mean1, sd1, rng_base);
        data[i+half] = stan::math::normal_rng(mean2, sd2, rng_base);
    }

    HypersFixedNNIG hy(4, 1.5, 2.0, 2.0);
    SimpleMixture mix(1.0);
    //Neal2<NNIGHierarchy, HypersFixedNNIG, SimpleMixture> sampler2(
    //    data, mix, hy);
    Neal8<NNIGHierarchy, HypersFixedNNIG, SimpleMixture> sampler8(
        data, 3, mix, hy);

    // Run samplers
    //sampler2.run();
    sampler8.run();

    // Density stuff
    std::vector<double> grid = {1,1.5,2,2.5,3,3.5,4};
    sampler8.eval_density(grid);
    //sampler8.write_density_to_file();

    //sampler2.eval_density(grid);
    //sampler2.write_density_to_file();

    // Clustering stuff
    //unsigned int i_cap = sampler8.cluster_estimate();
    //std::cout << "Best clustering: at iteration " << i_cap << std::endl;
    //sampler8.write_final_clustering_to_file();
    //sampler8.write_best_clustering_to_file();
    //sampler8.write_chain_to_file();

    return 0;
}
