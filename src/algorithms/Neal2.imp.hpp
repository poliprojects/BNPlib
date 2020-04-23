#ifndef NEAL2_IMP_HPP
#define NEAL2_IMP_HPP

#include "Neal2.hpp"

template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Neal2<Hierarchy, Hypers, Mixture>::print_startup_message(){
    std::cout << "Running Neal2 algorithm..." << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy, Hypers, Mixture>::initialize(){
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, this->num_clusters);

    for(int h = 0; h < this->num_clusters; h++){
      this->allocations.push_back(h);
    }
    for(int j = this->num_clusters; j < this->data.rows(); j++){
        int num = distribution(generator);
        this->allocations[j] = num;
    }
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy, Hypers, Mixture>::sample_allocations(){
    unsigned int k, n_unique, singleton;
    unsigned int n = this->data.rows();

    for(int i = 0; i < n; i++){ // for each data unit data[i]

        // Initialize cardinalities of unique values
        std::vector<int> card(this->unique_values.size(), 0);      
        for(int j = 0; j < n; j++){
            card[ this->allocations[j] ] += 1;
        }

        singleton = 0;
        n_unique = this->unique_values.size();

        if(card[ this->allocations[i] ] == 1){ // datum i is a singleton
            singleton = 1;
        }
        
        card[ this->allocations[i] ] -= 1;

        
        // Draw a NEW value for ci
        Eigen::VectorXd probas(n_unique+(1-singleton)); 
        

        auto M = this->mixture.get_totalmass();
        double tot = 0.0;
        
        for(int k = 0; k < n_unique; k++){
            probas(k) = this->mixture.prob_existing_cluster(card[k],n) *
            	this->unique_values[k].like(this->data.row(i))(0);

            if(singleton == 1 && k == i){
              probas(i) = this->mixture.prob_new_cluster(n, n_unique) *
                this->unique_values[0].eval_marg(this->data.row(i))(0);
            }
 
            tot += probas(k);
        }

        if(singleton == 0){
            probas(n_unique) = this->mixture.prob_new_cluster(n, n_unique) *
                this->unique_values[0].eval_marg(this->data.row(i))(0);
            tot += probas(n_unique);
        }

        // Normalize
        probas = probas / tot;
        
        unsigned int c_new = stan::math::categorical_rng(probas, this->rng) - 1;
        
        if(singleton == 1){
            if(c_new == this->allocations[i]){
                // case 1 of 4: SINGLETON - SINGLETON
                Eigen::VectorXd temp;
                temp=this->data.row(i); // initialize with datum if univariate,
                                        // with a vector (dim p) if multi
                this->unique_values[ this->allocations[i]
                    ].sample_given_data(temp);
                
            }
            else{ // case 2 of 4: SINGLETON - CLUSTER
                this->unique_values.erase(
                    this->unique_values.begin() + this->allocations[i]);
                
                int tmp = this->allocations[i];
                this->allocations[i] = c_new;
                for(auto &c : this->allocations){ // relabeling
                    if(c > tmp){
                        c -= 1;
                    }
                }
            } // end of else
        } // end of if(singleton == 1)

        else{ // if singleton == 0
            if(c_new == n_unique){ // case 3 of 4: NOT SINGLETON - SINGLETON

                Hierarchy<Hypers> new_unique(
                    this->unique_values[0].get_hypers());
		Eigen::VectorXd temp;
                temp=this->data.row(i);
                new_unique.sample_given_data(temp);
                this->unique_values.push_back(new_unique); 
                this->allocations[i] = n_unique;
            }
            else{ // case 4 of 4: NOT SINGLETON - CLUSTER
                this->allocations[i] = c_new;
            }

        } // end of else

    } // end of for(int i = 0; i < n; i++) loop

} // end of sample_allocations()


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy, Hypers, Mixture>::sample_unique_values(){

    this->num_clusters = this->unique_values.size();
    std::vector<std::vector<unsigned int>> clust_idxs(this->num_clusters);
    unsigned int n = this->allocations.size();
    for(int i = 0; i < n; i++){ // save different cluster in each row
        clust_idxs[ this->allocations[i] ].push_back(i);
    }

    for(int j = 0; j < this->num_clusters; j++){
	Eigen::MatrixXd curr_data(this->data.rows(), this->data.cols());
     	int k=0;

        for(auto &idx : clust_idxs[j]){
            curr_data.row(k)=this->data.row(idx);	
            k+=1;
	}
        curr_data.conservativeResize(k, Eigen::NoChange);
            // TODO: più efficiente?
        this->unique_values[j].sample_given_data(curr_data);
    }

}



template<template <class> class Hierarchy, class Hypers, class Mixture>
Eigen::VectorXd Neal2<Hierarchy, Hypers, Mixture>::eval_density_specific(Hierarchy<Hypers> &temp_hier,unsigned int n){
    double M = this->mixture.get_totalmass();
    Eigen::VectorXd dens_addendum(this->density.first.rows());

        // Component from G0 (exploit conjugacy using explicit expression)
    dens_addendum = M * temp_hier.eval_marg(this->density.first) / (M+n); 
    
    return dens_addendum;
}


#endif // NEAL2_IMP_HPP
