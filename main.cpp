#include <iostream>
#include  <fstream>

//#include <boost/random/random_number_generator.hpp>
//#include <boost/random/detail/qrng_base.hpp>

#include "includes_main.hpp"
#include "math.h"

int main(){
    unsigned int n = 100;
    std::vector<double> data(n);
    unsigned int half = data.size()/2;

    double mean1 = 4.0;
    double mean2 = 7.0;
    double sd1   = 1.0;
    double sd2   = 1.0;

    std::default_random_engine generator;
    std::normal_distribution<double> N1(mean1,sd1);
    std::normal_distribution<double> N2(mean2,sd2);

    for(int i = 0; i < half; i++){
        data[i]      = N1(generator);
        data[i+half] = N2(generator);
    }

    //std::ofstream file;
    //file.open("data.csv");
    //for(auto &d : data){
    //    file << d << ",";
    //}
    //file << std::endl;
    //file.close();

  
    HypersFixedNNIG hy(5.0, 1.0, 2.0, 2.0); // mu0, lambda, alpha0, beta0

   
    SimpleMixture mix(0.25); // total mass
    Neal2<HierarchyNNIG, HypersFixedNNIG, SimpleMixture> sampler2(
        data, mix, hy);
    //Neal8<HierarchyNNIG, HypersFixedNNIG, SimpleMixture> sampler8(
    //    data, 3, mix, hy);
	
    // Run samplers
    sampler2.run();
    //sampler8.run();

    
	



    // Density stuff
    std::vector<double> grid;
    double temp = 0.0;
    double step = 0.05;
    double upp_bnd = 10.0;
    while(temp <= upp_bnd){
        grid.push_back(temp);
        temp += step;
    }
    //sampler8.eval_density(grid);
    //sampler8.write_density_to_file("densityfinal0.25.csv");

    sampler2.eval_density(grid);
    sampler2.write_density_to_file("densityneal2.csv");
	//unsigned int i_cap = sampler2.cluster_estimate();
    //std::cout << "Best clustering: at iteration " << i_cap << std::endl;
    //sampler2.write_final_clustering_to_file();
    //sampler2.write_best_clustering_to_file();

    // Clustering stuff
    //unsigned int i_cap = sampler8.cluster_estimate();
    //std::cout << "Best clustering: at iteration " << i_cap << std::endl;
    //sampler8.write_final_clustering_to_file("clust_final0.25.csv");
    //sampler8.write_best_clustering_to_file("clust_best0.25.csv");
    //sampler8.write_chain_to_file();

    return 0;
}
