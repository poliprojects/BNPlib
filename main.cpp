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

	std::vector<double> v;
	
	while(std::getline(iss, str2, ',')){
		double val = ::atof(str2.c_str());
		v.push_back(val);
		
	}
	file.close();
    Eigen::VectorXd data = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(v.data(), v.size());


    double mu0, lambda, alpha0, beta0;
    //std::cout << "Insert mu0, lambda, alpha0, beta0 values:" << std::endl;
    //std::cin >> mu0 >> lambda >> alpha0 >> beta0; // 5.0 0.1 2.0 2.0
    // HypersFixedNNIG hy(mu0, lambda, alpha0, beta0);

    double totalmass;
    //std::cout << "Insert total mass value:" << std::endl; 
    //std::cin >> totalmass; //1.0
    // DirichletMixture mix(totalmass);

    int n_aux(3);
    //std::cout << "Insert number of auxiliary blocks:" << std::endl;
    //std::cin >> n_aux;

    //std::ofstream file;
    //file.open("data.csv");
    //for(auto &d : data){
    //    file << d << ",";
    //}
    //file << std::endl;
    //file.close();
 

    HypersFixedNNIG hy(5.0, 1.0, 2.0, 2.0); // mu0, lambda, alpha0, beta0
    DirichletMixture mix(1); // total mass

    Neal8<HierarchyNNIG, HypersFixedNNIG, DirichletMixture> sampler(
        data, n_aux, mix, hy);
	
    // Run sampler
    sampler.run();

    return 0;

    // Density stuff
    Eigen::VectorXd grid;
    double temp = 0.0;
    double step = 0.05;
    double upp_bnd = 10.0;
    while(temp <= upp_bnd){
        grid<<temp;
        temp += step;
    }
    //sampler.eval_density(grid);
    //sampler.write_density_to_file("density_m50.csv");

    // Density and clustering stuff
    //sampler.eval_density(grid);
    //sampler.write_density_to_file();
	//unsigned int i_cap = sampler.cluster_estimate();
    //std::cout << "Best clustering: at iteration " << i_cap << std::endl;
    //sampler.write_final_clustering_to_file();
    //sampler.write_best_clustering_to_file();

    return 0;
}
