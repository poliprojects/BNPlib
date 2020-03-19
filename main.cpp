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

	Eigen::VectorXd data;
	while(std::getline(iss, str2, ',')){
		double val = ::atof(str2.c_str());
		data<<val;
	}
	file.close();

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
<<<<<<< HEAD
    //Neal2<HierarchyNNIG, HypersFixedNNIG, DirichletMixture> sampler2(
      //  data, mix, hy);
    Neal8<HierarchyNNIG, HypersFixedNNIG, DirichletMixture> sampler8(
        data, n_aux, mix, hy);
	
    // Run sampler(s)
    //sampler2.run();
    //sampler8.run();
=======
    Neal8<HierarchyNNIG, HypersFixedNNIG, DirichletMixture> sampler(
        data, mix, hy);
    //Neal8<HierarchyNNIG, HypersFixedNNIG, DirichletMixture> sampler(
      //  data, 3, mix, hy);
	
    // Run sampler(s)
    sampler.run();
>>>>>>> 3a5ab117e3f3e9471e7c90d1f503ec865c5d9b81

    //Neal8<HierarchyNNIG, HypersFixedNNIG, DirichletMixture> sampler(
      // 	data, n_aux, mix, hy);
	
    //sampler.run();


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

<<<<<<< HEAD
    //sampler2.eval_density(grid);
    //sampler2.write_density_to_file("densityneal2.csv");
	//unsigned int i_cap = sampler2.cluster_estimate();
=======
    sampler.eval_density(grid);
    //sampler.write_density_to_file("densityneal2.csv");
	//unsigned int i_cap = sampler.cluster_estimate();
>>>>>>> 3a5ab117e3f3e9471e7c90d1f503ec865c5d9b81
    //std::cout << "Best clustering: at iteration " << i_cap << std::endl;
    //sampler.write_final_clustering_to_file();
    //sampler.write_best_clustering_to_file();

    // Clustering stuff
    //unsigned int i_cap = sampler.cluster_estimate();
    //std::cout << "Best clustering: at iteration " << i_cap << std::endl;
    //sampler.write_final_clustering_to_file("clust_final0.25.csv");
    //sampler.write_best_clustering_to_file("clust_best1.csv");
    //sampler.write_chain_to_file();

    //return 0;
}
