#ifndef ALGORITHM_IMP_HPP
#define ALGORITHM_IMP_HPP

#include "Algorithm.hpp"


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Algorithm<Hierarchy, Hypers, Mixture>::save_iteration(unsigned int iter){
    //std::cout << "Iteration # " << iter << " / " << maxiter-1 <<
    //    std::endl; // DEBUG
    IterationOutput iter_out;
    *iter_out.mutable_allocations() = {allocations.begin(), allocations.end()};

    for(int i = 0; i < unique_values.size(); i++){
        UniqueValues temp;
        for(auto &par : unique_values[i].get_state()){
            temp.add_params(par);
        }
        iter_out.add_uniquevalues();
        *iter_out.mutable_uniquevalues(i) = temp;
    }

    chain.add_chain();
    *chain.mutable_chain(iter-burnin) = iter_out;

    //print_state(); //DEBUG
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Algorithm<Hierarchy, Hypers, Mixture>::print_state(){
    for (int h = 0; h < num_clusters; h++) {
        std::cout << "Parameters: ";

        for (auto c : unique_values[h].get_state()){
            std::cout << c << " " << std::endl;
        }
        std::cout << std::endl;
    }
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Algorithm<Hierarchy, Hypers, Mixture>::print_ending_message(){
    std::cout << "Done" << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
unsigned int Algorithm<Hierarchy, Hypers, Mixture>::cluster_estimate(){
    // also returns the index of the estimate in the chain object

    unsigned int niter = maxiter - burnin;
    Eigen::VectorXd errors(niter);
    int n = data.size();
    Eigen::MatrixXd tot_diss(n, n);
    tot_diss = Eigen::MatrixXd::Zero(n, n);
    std::vector<Eigen::MatrixXd> all_diss;
    IterationOutput temp;
    
    for(int h = 0; h < niter; h++){
        temp = *chain.mutable_chain(h);
        Eigen::MatrixXd dissim(n, n);
        dissim = Eigen::MatrixXd::Zero(n, n);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < i; j++){
                if(temp.allocations(i) == temp.allocations(j)){
                    dissim(i,j) = 1;
                }
            }
        }
    all_diss.push_back(dissim);
    tot_diss = tot_diss + dissim;
    }

    tot_diss = tot_diss / niter;

    for(int h = 0; h < niter; h++){
        // Compute error in Frobenius norm
        errors(h) = (tot_diss-all_diss[h]).norm();
    }
    
    //std::cout << errors << std::endl; // DEBUG
    std::ptrdiff_t i;
    int min_err = errors.minCoeff(&i);

    //std::ofstream file;
    //file.open("dissim_matr_mean.csv");
    //file << tot_diss;
    //file.close();

    //std::ofstream file2;
    //file2.open("dissim_matr_best.csv");
    //file2 << all_diss[i];
    //file2.close();

    best_clust = chain.chain(i);
    std::cout << best_clust.uniquevalues_size() <<
        " clusters were found via least square minimization" << std::endl;
    return i;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Algorithm<Hierarchy, Hypers, Mixture>::write_final_clustering_to_file(
        std::string filename){
    // number,datum,cluster,params1,params2,...
    
    std::ofstream file;
    file.open(filename);

    for(int i = 0; i < data.size(); i++){
        auto params = unique_values[ allocations[i] ].get_state();
        file << i << "," << data[i] << "," << allocations[i];
        for(int j = 0; j < params.size(); j++){
            file << "," << params[j];
        }
        file << std::endl;
    }
    file.close();
    std::cout << "Succesfully wrote to " << filename << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Algorithm<Hierarchy, Hypers, Mixture>::write_best_clustering_to_file(
    std::string filename){
    // number,datum,cluster,params1,params2,...

    std::ofstream file;
    file.open(filename);

    for(int i = 0; i < data.size(); i++){
        unsigned int ci = best_clust.allocations(i);
        file << i << "," << data[i] << "," << ci;
        for(int j = 0; j < best_clust.uniquevalues(ci).params_size(); j++){
            file << "," << best_clust.uniquevalues(ci).params(j);
        }
        file << std::endl;
    }
    file.close();
    std::cout << "Succesfully wrote to " << filename << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Algorithm<Hierarchy, Hypers, Mixture>::write_chain_to_file(
    std::string filename){
    // number,datum,cluster,params1,params2,...

    std::ofstream file;
    file.open(filename);

    // for each iteration of the algorithm
    for(int iter = 0; iter < chain.chain_size(); iter++){
        // for each data point
        for(int i = 0; i < data.size(); i++){
            auto state_iter = chain.chain(iter);
            unsigned int ci = state_iter.allocations(i);
            file << iter << "," << i << "," << data[i] << "," << ci;
            for(int j = 0; j < state_iter.uniquevalues(ci).params_size(); j++){
                file << "," << state_iter.uniquevalues(ci).params(j);
            }
            file << std::endl;
        }
    }

    file.close();
    std::cout << "Succesfully wrote to " << filename << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Algorithm<Hierarchy, Hypers, Mixture>::write_density_to_file(
    std::string filename){
    std::ofstream file;
    file.open(filename);

    for(int i = 0; i < density.first.size(); i++){
        file << density.first[i] << "," << density.second(i) << std::endl;
    }
    
    file.close();
    std::cout << "Succesfully wrote to " << filename << std::endl;
}


#endif // ALGORITHM_IMP_HPP
