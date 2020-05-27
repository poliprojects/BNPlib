#ifndef NEAL2_IMP_HPP
#define NEAL2_IMP_HPP

#include "Neal2.hpp"

// Algorithm functions
template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Neal2<Hierarchy, Hypers, Mixture>::print_startup_message(){
    std::cout << "Running Neal2 algorithm..." << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy, Hypers, Mixture>::initialize(){
    cardinalities.reserve(data.rows());
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distro(0,this->init_num_clusters-1);

    for(size_t h = 0; h < this->init_num_clusters; h++){
      allocations.push_back(h);
      cardinalities.push_back(1);
    }
    for(size_t j = this->init_num_clusters; j < data.rows(); j++){
        unsigned int clust = distro(generator);
        allocations[j] = clust;
        cardinalities[clust] += 1;
    }
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy, Hypers, Mixture>::sample_allocations(){
    unsigned int n_unique, singleton;
    unsigned int n = data.rows();
    double M = this->mixture.get_totalmass();

    for(size_t i = 0; i < n; i++){
        // for each data unit data[i]
        Eigen::Matrix<double, 1, Eigen::Dynamic> datum = data.row(i);

        // Initialize cardinalities of unique values
        singleton = 0;
        n_unique = unique_values.size();

        if(cardinalities[ allocations[i] ] == 1){
            singleton = 1; // datum i is a singleton
        }
        
        // Remove datum from cluster
        cardinalities[ allocations[i] ] -= 1;
        
        // Compute probabilities of clusters
        Eigen::VectorXd probas(n_unique+(1-singleton));

        double tot = 0.0;
        
        for(size_t k = 0; k < n_unique; k++){
            probas(k) = this->mixture.prob_existing_cluster(
                cardinalities[k], n) * unique_values[k].like(datum)(0);

            if(singleton == 1 && k == i){
                probas(i) = this->mixture.prob_new_cluster(n, n_unique) *
                    unique_values[0].eval_marg(datum)(0);
            }
 
            tot += probas(k);
        }

        if(singleton == 0){
            probas(n_unique) = this->mixture.prob_new_cluster(n, n_unique) *
                unique_values[0].eval_marg(datum)(0);
            tot += probas(n_unique);
        }



        // Normalize
        probas = probas / tot;
        
        // Draw a NEW value for ci
        unsigned int c_new = stan::math::categorical_rng(probas, this->rng) - 1;
        
        // Assign datum to the new cluster
        if(singleton == 1){
            if(c_new == allocations[i]){
                // case 1 of 4: SINGLETON - SINGLETON

                unique_values[ allocations[i] ].sample_given_data(datum);
                cardinalities[c_new] += 1;
            }

            else{ // case 2 of 4: SINGLETON - CLUSTER

                unique_values.erase( unique_values.begin() + allocations[i] );
                
                unsigned int c_old = allocations[i];
                allocations[i] = c_new;
                for(auto &c : allocations){ // relabeling
                    if(c > c_old){
                        c -= 1;
                    }
                }
                cardinalities[c_new] += 1;
                cardinalities.erase(cardinalities.begin() + c_old);
            } // end of else
        } // end of if(singleton == 1)

        else{ // if singleton == 0
            if(c_new == n_unique){ // case 3 of 4: NOT SINGLETON - SINGLETON

                Hierarchy<Hypers> new_unique( unique_values[0].get_hypers() );

                new_unique.sample_given_data(datum);
                unique_values.push_back(new_unique);
                allocations[i] = n_unique;
                cardinalities.push_back(1);
            }
            else{ // case 4 of 4: NOT SINGLETON - CLUSTER

                allocations[i] = c_new;
                cardinalities[c_new] += 1;
            }
        } // end of else
    } // end of datum[i] loop
} // end of sample_allocations()


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy, Hypers, Mixture>::sample_unique_values(){
    unsigned int n_unique = unique_values.size();

    std::vector<std::vector<unsigned int>> clust_idxs(n_unique);
    unsigned int n = allocations.size();
    for(size_t i = 0; i < n; i++){ // save different cluster in each row
        clust_idxs[ allocations[i] ].push_back(i);
    }

    for(size_t i = 0; i < n_unique; i++){
        unsigned int curr_size = clust_idxs[i].size();
        Eigen::MatrixXd curr_data(curr_size, data.cols());
        for(size_t j = 0; j < curr_size; j++){
            curr_data.row(j) = data.row( clust_idxs[i][j] );
        }
        unique_values[i].sample_given_data(curr_data);
    }
}


// Other tools
template<template <class> class Hierarchy, class Hypers, class Mixture>
Eigen::VectorXd Neal2<Hierarchy, Hypers, Mixture>::density_marginal_component(
    Hierarchy<Hypers> &temp_hier, unsigned int n){
    // Component from G0 (exploit conjugacy using explicit expression)
    double M = this->mixture.get_totalmass();
    Eigen::VectorXd dens_addendum(this->density.first.rows());
    dens_addendum = temp_hier.eval_marg(this->density.first) * M/(M+n);
    return dens_addendum;
}


#endif // NEAL2_IMP_HPP
