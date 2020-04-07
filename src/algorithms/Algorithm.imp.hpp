#ifndef ALGORITHM_IMP_HPP
#define ALGORITHM_IMP_HPP

#include "Algorithm.hpp"

template<template <class> class Hierarchy, class Hypers, class Mixture>
Eigen::MatrixXd Algorithm<Hierarchy, Hypers, Mixture>::proto_param_to_matrix(const Param &par){
    Eigen::MatrixXd Par_matrix= Eigen::MatrixXd::Zero(par.par_cols_size(), par.par_cols(0).elems_size());
    for(int h = 0; h < par.par_cols_size(); h++){
        for(int j = 0; j < par.par_cols(h).elems_size(); j++){
        Par_matrix(j,h)=par.par_cols(h).elems(j);
        }
    }

return Par_matrix;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
IterationOutput Algorithm<Hierarchy, Hypers, Mixture>::get_state_as_proto(unsigned int iter){

    IterationOutput iter_out;
    *iter_out.mutable_allocations() = {allocations.begin(), allocations.end()};

    for(int i = 0; i < unique_values.size(); i++){
        UniqueValues Uniquevalues_temp;
        for(int k = 0; k< unique_values[i].get_state().size(); k++){
            Eigen::MatrixXd par_temp=unique_values[i].get_state()[k];
            Param par_temp_proto;
            for(int j=0; j<par_temp.cols(); j++){
                Par_Col col_temp;
                for(int h=0; h<par_temp.rows(); h++){
                col_temp.add_elems(par_temp(h,j));
                }
                par_temp_proto.add_par_cols();
                *par_temp_proto.mutable_par_cols(j)=col_temp;     
            }
            Uniquevalues_temp.add_params();
            *Uniquevalues_temp.mutable_params(k)=par_temp_proto;
        }
        iter_out.add_uniquevalues();
        *iter_out.mutable_uniquevalues(i) = Uniquevalues_temp;
    }
    return iter_out;
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
unsigned int Algorithm<Hierarchy, Hypers, Mixture>::cluster_estimate(BaseCollector* collector){
    // also returns the index of the estimate in the chain object

    unsigned int niter = maxiter - burnin;
    Eigen::VectorXd errors(niter);
    int n = data.rows();
    Eigen::MatrixXd tot_diss(n, n);
    tot_diss = Eigen::MatrixXd::Zero(n, n);
    std::vector<Eigen::MatrixXd> all_diss;
    IterationOutput temp;
    
    for(int h = 0; h < niter; h++){
        temp = collector->get_chains()[h];
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

    best_clust = collector->get_chains()[i];
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

    for(int i = 0; i < data.rows(); i++){ // TODO tutti i write da modificare con dati multivar
        auto params = unique_values[ allocations[i] ].get_state();
        file << i << "," << data.row(i) << "," << allocations[i];
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

    for(int i = 0; i < data.rows(); i++){
        unsigned int ci = best_clust.allocations(i);
        file << i << "," << data.row(i) << "," << ci;
        for(int j = 0; j < best_clust.uniquevalues(ci).params_size(); j++){
            file << "," << this->proto_param_to_matrix(best_clust.uniquevalues(ci).params(j));
        }
        file << std::endl;
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
