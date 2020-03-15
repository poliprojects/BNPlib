#include <iostream>
#include <fstream>

//#include <boost/random/random_number_generator.hpp>
//#include <boost/random/detail/qrng_base.hpp>

#include "includes_main.hpp"
#include "math.h"

int main(int argc, char *argv[]){
	// Read data from main arg
	std::ifstream file;
	file.open(argv[1]);
	std::string str, str2;
	std::getline(file, str);
	std::istringstream iss(str);

	std::vector<double> data;
	while(std::getline(iss, str2, ',')){
		double val = ::atof(str2.c_str());
		data.push_back(val);
	}
	file.close();

    double mu0, lambda, alpha0, beta0;
    std::cout << "Insert mu0, lambda, alpha0, beta0 values:" << std::endl;
    std::cin >> mu0 >> lambda >> alpha0 >> beta0; // 5.0 0.1 2.0 2.0
    HypersFixedNNIG hy(mu0, lambda, alpha0, beta0);

    double totalmass;
    std::cout << "Insert total mass value:" << std::endl; 
    std::cin >> totalmass; //1.0
    DirichletMixture mix(totalmass);

    unsigned int n_aux;
    std::cout << "Insert number of auxiliary blocks:" << std::endl;
    std::cin >> n_aux;

    Neal8<HierarchyNNIG, HypersFixedNNIG, DirichletMixture> sampler(
       	data, n_aux, mix, hy);
	
    sampler.run();

    return 0;

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
    //sampler8.write_density_to_file("density_m50.csv");

    //sampler2.eval_density(grid);
    //sampler2.write_density_to_file("densityneal2.csv");
	//unsigned int i_cap = sampler2.cluster_estimate();
    //std::cout << "Best clustering: at iteration " << i_cap << std::endl;
    //sampler2.write_final_clustering_to_file();
    //sampler2.write_best_clustering_to_file();

    // Clustering stuff
    //unsigned int i_cap = sampler8.cluster_estimate();
    //std::cout << "Best clustering: at iteration " << i_cap << std::endl;
    //sampler8.write_final_clustering_to_file("clust_final0.25.csv");
    //sampler8.write_best_clustering_to_file("clust_best1.csv");
    //sampler8.write_chain_to_file();

    //return 0;
}
