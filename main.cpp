#include <iostream>
#include <fstream>

//#include <boost/random/random_number_generator.hpp>
//#include <boost/random/detail/qrng_base.hpp>
#include <chrono>

#include "includes_main.hpp"
#include "math.h"

int main(int argc, char *argv[]){
// 3D-vectorial data
    Eigen::MatrixXd data(3,5);
    data << 1.3, 0.9, 8.8, 2.0, -1.3,
            2.3, 5.1, 4.4, 0.0, -3.2,
            3.0, 4,0, 3.3, 1.1, +1.5;

    Eigen::VectorXd mu0(3);
    mu0 << 3.0, 3.0, 3.0;
    Eigen::MatrixXd lambda0 = 2 * Eigen::Matrix<double, 3, 3>::Identity();
    double totalmass = 1.0;
    int n_aux = 3;
    HypersDummy hy(mu0, lambda0);
    DirichletMixture mix(totalmass); // total mass
    Neal8<HierarchyDummy, HypersDummy, DirichletMixture> sampler(
        data, n_aux, mix, hy);

    BaseCollector *f;
    std::string collector(argv[2]);
    if(collector=="FileCollector"){
        std::string filename(argv[3]);
        f=new FileCollector(filename);}
    if(collector=="MemoryCollector"){
        f=new MemoryCollector();}
   
    sampler.run(f);

    return 0;





//    // PROVA TIME //////////////////////
//    
//    Eigen::MatrixXd data_test=Eigen::MatrixXd::Identity(1000,2);
//    // try with diff cluster size
//    //std::vector<int> clust_idxs={0,30,20,5,6,7,22,34,91,4,35,76,90};
//    //std::vector<int> clust_idxs={1};
//    std::vector<int> clust_idxs(1000) ; 
//
//    std::iota (std::begin(clust_idxs), std::end(clust_idxs), 0); // Fill with 0, 1, ..., 999.
//
//    std::chrono::time_point<std::chrono::system_clock> start_1, end_1,start_2, end_2;
//
// 
//
//    // FIRST WAY 
//    start_1 = std::chrono::system_clock::now();    
//    Eigen::MatrixXd curr_data(data_test.rows(), data_test.cols());
//     	int k=0;
//
//
//        for(auto &idx : clust_idxs){
//            curr_data.row(k)=data_test.row(idx);	
//            k+=1;
//	}
//        curr_data.conservativeResize(k,Eigen::NoChange); 
//    end_1 = std::chrono::system_clock::now();
//        
//    // SECOND WAY
//    start_2 = std::chrono::system_clock::now();    
//    std::vector<double> curr_data_2;
//
//        for(auto &idx : clust_idxs){
//            for (int i=0; i<data_test.cols(); i++){
//                curr_data_2.push_back(data_test.row(idx)(i));
//            }
//	}
//	
//
//            
//    Eigen::MatrixXd curr_data_map = Eigen::Map<Eigen::MatrixXd>(curr_data_2.data(),curr_data_2.size()/data_test.cols() ,data_test.cols()); 
//    end_2 = std::chrono::system_clock::now();
//    typedef std::chrono::duration<int, std::ratio<1, 100000000>> shakes;
//
//    int elapsed_seconds_1 = std::chrono::duration_cast<shakes>(end_1-start_1).count();
//    int elapsed_seconds_2 = std::chrono::duration_cast<shakes>(end_2-start_2).count();
//    std::cout<<"first time :"<<elapsed_seconds_1<<std::endl;
//    std::cout<<"second time :"<<elapsed_seconds_2<<std::endl;
//	return 0;
//	// END PROVA TIME /////////////////////////
//	
//	// Read data from main arg
//	std::ifstream file;
//	file.open(argv[1]);
//	std::string str, str2;
//	std::getline(file, str);
//	std::istringstream iss(str);
//
//	std::vector<double> v;
//	
//	while(std::getline(iss, str2, ',')){
//		double val = ::atof(str2.c_str());
//		v.push_back(val);
//		
//	}
//	file.close();
//    Eigen::VectorXd data = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(v.data(), v.size()); // TODO: meglio con conservative resize?
//
//
//    double mu0, lambda, alpha0, beta0;
//    //std::cout << "Insert mu0, lambda, alpha0, beta0 values:" << std::endl;
//    //std::cin >> mu0 >> lambda >> alpha0 >> beta0; // 5.0 0.1 2.0 2.0
//    // HypersFixedNNIG hy(mu0, lambda, alpha0, beta0);
//
//    double totalmass;
//    //std::cout << "Insert total mass value:" << std::endl; 
//    //std::cin >> totalmass; //1.0
//    // DirichletMixture mix(totalmass);
//
//    int n_aux(3);
//    //std::cout << "Insert number of auxiliary blocks:" << std::endl;
//    //std::cin >> n_aux;
//
//    //std::ofstream file;
//    //file.open("data.csv");
//    //for(auto &d : data){
//    //    file << d << ",";
//    //}
//    //file << std::endl;
//    //file.close();
// 
//
//    HypersFixedNNIG hy(5.0, 1.0, 2.0, 2.0); // mu0, lambda, alpha0, beta0
//    DirichletMixture mix(1); // total mass
//
//    //Neal8<HierarchyNNIG, HypersFixedNNIG, DirichletMixture> sampler(
//      //  data, n_aux, mix, hy);
//    Neal2<HierarchyNNIG, HypersFixedNNIG, DirichletMixture> sampler2(
//        data, mix, hy);
//	
//    // Run sampler
//	BaseCollector *f;
//    std::string collector(argv[2]);
//    
//    if(collector=="FileCollector"){
//        std::string filename(argv[3]);
//        f=new FileCollector(filename);}
//    if(collector=="MemoryCollector"){
//        f=new MemoryCollector();}
//   
//    
//    sampler2.run(f);
//    //sampler.run(f);
//
//
//    // Density stuff
//    //Eigen::VectorXd grid;
//    double temp = 0.0;
//    double step = 0.05;
//    double upp_bnd = 10.0;
//    std::vector<double> v_temp;
//    while(temp <= upp_bnd){
//        v_temp.push_back(temp);
//        temp += step;
//    }
//    Eigen::VectorXd grid = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(v_temp.data(), v_temp.size()); 
// 
// 
//   
//    //sampler.eval_density(grid, f);
//    //sampler.write_density_to_file("density_m50.csv");
//	//unsigned int i_cap = sampler.cluster_estimate(f);
//    //std::cout << "Best clustering: at iteration " << i_cap << std::endl;
//    //sampler.write_final_clustering_to_file();
//    //sampler.write_best_clustering_to_file();
//
//
//    //sampler2.eval_density(grid,f);
//    //sampler2.write_density_to_file();
//	//unsigned int i_cap = sampler2.cluster_estimate(f);
//    //std::cout << "Best clustering: at iteration " << i_cap << std::endl;
//    //sampler2.write_final_clustering_to_file();
//    //sampler2.write_best_clustering_to_file();
//
//    return 0;
}
