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
    std::uniform_int_distribution<int> distribution(0,num_clusters);

    for(int h = 0; h < num_clusters; h++){
      allocations.push_back(h);
    }
    for(int j = num_clusters; j < data.size(); j++){
        int num = distribution(generator); //TODO da stan?
        allocations[j] = num;
    }
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy, Hypers, Mixture>::sample_allocations(){
    unsigned int k, n_unique, singleton;
    unsigned int n = data.size();

    for(int i = 0; i < n; i++){ // for each data unit data[i]

        // Initialize cardinalities of unique values
        std::vector<int> card(unique_values.size(), 0);      
        for(int j = 0; j < n; j++){
            card[ allocations[j] ] += 1;
        }

        singleton = 0;
        n_unique = unique_values.size();

        if(card[ allocations[i] ] == 1){ // datum i is a singleton
            singleton = 1;
        }
        
        card[ allocations[i] ] -= 1;

        
        // Draw a NEW value for ci
        Eigen::VectorXd probas(n_unique+(1-singleton)); 
        

        auto M = mixture.get_totalmass();
        double tot = 0.0;

        for(int k = 0; k < n_unique; k++){
            probas(k) = card[k] * unique_values[k].like(data[i]) / (
                n-1+M);
            if(singleton == 1 && k == i){
                std::shared_ptr<Hypers> hy = unique_values[0].get_hypers();
                double mu0    = hy->get_mu0();
                double lambda = hy->get_lambda();
                double alpha0 = hy->get_alpha0();
                double beta0  = hy->get_beta0();
                
                double sigtilde = sqrt( beta0*(lambda+1)/(alpha0*lambda) );
                probas(i,0) = M * exp(stan::math::student_t_lpdf(data[i],
                    2*alpha0, mu0, sigtilde)) / (n-1+M);
            } 
            tot += probas(k);
        }

        if(singleton == 0){
            std::shared_ptr<Hypers> hy = unique_values[0].get_hypers();
                double mu0    = hy->get_mu0();
                double lambda = hy->get_lambda();
                double alpha0 = hy->get_alpha0();
                double beta0  = hy->get_beta0();

            double sigtilde = sqrt( beta0*(lambda+1)/(alpha0*lambda) );
            probas(n_unique,0) = M * exp(stan::math::student_t_lpdf(data[i],
                2*alpha0, mu0, sigtilde)) / (n-1+M);
            tot += probas(n_unique,0);
        }

        // Normalize
        probas = probas / tot;
        
        unsigned int c_new = stan::math::categorical_rng(probas, rng) - 1;
        
        if(singleton == 1){
            if(c_new == allocations[i]){ // case 1 of 4: SINGLETON - SINGLETON
                std::shared_ptr<Hypers> hy = unique_values[0].get_hypers();
                double mu0    = hy->get_mu0();
                double lambda = hy->get_lambda();
                double alpha0 = hy->get_alpha0();
                double beta0  = hy->get_beta0();    

                std::vector<double> par_pair(2);
                double sigma2_new = stan::math::inv_gamma_rng(alpha0 + 1/2,
                    beta0+(lambda*pow(data[i]-mu0,2))/(2*(lambda+1)), rng);
                double mu_new = stan::math::normal_rng(
                    (lambda*mu0+data[i])/(lambda+1),
                    sqrt(sigma2_new/(lambda+1)), rng); 
                par_pair[0] = mu_new;
                par_pair[1] = sqrt(sigma2_new);
                unique_values[ allocations[i] ].set_state(par_pair);
                
            }
            else{ // case 2 of 4: SINGLETON - CLUSTER
                unique_values.erase( unique_values.begin()+allocations[i] );
                
                int tmp = allocations[i];
                allocations[i] = c_new;
                for(auto &c : allocations){ // relabeling
                    if(c > tmp){
                        c -= 1;
                    }
                }
            } // end of else
        } // end of if(singleton == 1)

        else{ // if singleton == 0
            if(c_new == n_unique){ // case 3 of 4: NOT SINGLETON - SINGLETON
                std::shared_ptr<Hypers> hy= unique_values[0].get_hypers();
                double mu0    = hy->get_mu0();
                double lambda = hy->get_lambda();
                double alpha0 = hy->get_alpha0();
                double beta0  = hy->get_beta0();

               std::vector<double> par_pair(2);

                double sigma2_new = stan::math::inv_gamma_rng(alpha0 + 1/2,
                    beta0+(lambda*pow(data[i]-mu0,2))/(2*(lambda+1)), rng);
                double mu_new = stan::math::normal_rng((lambda*mu0+data[i])/
                    (lambda+1), sqrt(sigma2_new/(lambda+1)), rng); 
                par_pair[0] = mu_new;
                par_pair[1] = sqrt(sigma2_new);
                Hierarchy<Hypers> new_unique(hy);
                new_unique.set_state(par_pair); 
                unique_values.push_back(new_unique); 
                allocations[i] = n_unique;
            }
            else{ // case 4 of 4: NOT SINGLETON - CLUSTER
                allocations[i] = c_new;
            }

        } // end of else

    } // end of for(int i = 0; i < n; i++) loop

} // end of sample_allocations()


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy, Hypers, Mixture>::sample_unique_values(){

    num_clusters = unique_values.size();
    std::vector<std::vector<unsigned int>> clust_idxs(num_clusters);
    unsigned int n = allocations.size();
    for(int i = 0; i < n; i++){ // save different cluster in each row
        clust_idxs[ allocations[i] ].push_back(i);
    }

    // DEBUG:
    //for(unsigned int j = 0; j < num_clusters; j++){ 
    //    std::cout << "Cluster #" << j << ": ";
    //    for(unsigned int i = 0; i < clust_idxs[j].size(); i++){
    //        std::cout << " " << clust_idxs[j][i];
    //    }
    //    std::cout << std::endl;
    //}

    for(int j = 0; j < num_clusters; j++){
        std::vector<double> curr_data;
        for(auto &idx : clust_idxs[j])
            curr_data.push_back(data[idx]);
        unique_values[j].sample_given_data(curr_data);
    }

    // std::cout << std::endl; // DEBUG
}

#endif // NEAL2_IMP_HPP
