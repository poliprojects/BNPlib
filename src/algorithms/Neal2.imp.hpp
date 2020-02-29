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
    for(int j = this->num_clusters; j < this->data.size(); j++){
        int num = distribution(generator); //TODO da stan?
        this->allocations[j] = num;
    }
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy, Hypers, Mixture>::sample_allocations(){
    unsigned int k, n_unique, singleton;
    unsigned int n = this->data.size();

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
            probas(k) = card[k] * this->unique_values[k].like(this->data[i]) / (
                n-1+M);
            if(singleton == 1 && k == i){
                std::shared_ptr<Hypers> hy=this->unique_values[0].get_hypers();
                double mu0    = hy->get_mu0();
                double lambda = hy->get_lambda();
                double alpha0 = hy->get_alpha0();
                double beta0  = hy->get_beta0();
                
                double sigtilde = sqrt( beta0*(lambda+1)/(alpha0*lambda) );
                probas(i,0) = M * exp(stan::math::student_t_lpdf(this->data[i],
                    2*alpha0, mu0, sigtilde)) / (n-1+M);
            } 
            tot += probas(k);
        }

        if(singleton == 0){
            std::shared_ptr<Hypers> hy = this->unique_values[0].get_hypers();
                double mu0    = hy->get_mu0();
                double lambda = hy->get_lambda();
                double alpha0 = hy->get_alpha0();
                double beta0  = hy->get_beta0();

            double sigtilde = sqrt( beta0*(lambda+1)/(alpha0*lambda) );
            probas(n_unique,0) = M * exp(stan::math::student_t_lpdf(
                this->data[i], 2*alpha0, mu0, sigtilde)) / (n-1+M);
            tot += probas(n_unique,0);
        }

        // Normalize
        probas = probas / tot;
        
        unsigned int c_new = stan::math::categorical_rng(probas, this->rng)-1;
        
        if(singleton == 1){
            if(c_new == this->allocations[i]){
                // case 1 of 4: SINGLETON - SINGLETON

                std::shared_ptr<Hypers> hy=this->unique_values[0].get_hypers();
                double mu0    = hy->get_mu0();
                double lambda = hy->get_lambda();
                double alpha0 = hy->get_alpha0();
                double beta0  = hy->get_beta0();    

                std::vector<double> par_pair(2);
                double sigma2_new = stan::math::inv_gamma_rng(alpha0 + 1/2,
                    beta0+(lambda*pow(this->data[i]-mu0,2))/(2*(lambda+1)),
                    this->rng);
                double mu_new = stan::math::normal_rng(
                    (lambda*mu0+this->data[i])/(lambda+1),
                    sqrt(sigma2_new/(lambda+1)), this->rng); 
                par_pair[0] = mu_new;
                par_pair[1] = sqrt(sigma2_new);
                this->unique_values[ this->allocations[i] ].set_state(par_pair);
                
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
                std::shared_ptr<Hypers> hy=this->unique_values[0].get_hypers();
                double mu0    = hy->get_mu0();
                double lambda = hy->get_lambda();
                double alpha0 = hy->get_alpha0();
                double beta0  = hy->get_beta0();

               std::vector<double> par_pair(2);

                double sigma2_new = stan::math::inv_gamma_rng(alpha0 + 1/2,
                    beta0+(lambda*pow(this->data[i]-mu0,2))/(2*(lambda+1)),
                    this->rng);
                double mu_new = stan::math::normal_rng((lambda*mu0+
                    this->data[i])/(lambda+1), sqrt(sigma2_new/(lambda+1)),
                    this->rng); 
                par_pair[0] = mu_new;
                par_pair[1] = sqrt(sigma2_new);
                Hierarchy<Hypers> new_unique(hy);
                new_unique.set_state(par_pair); 
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

    // DEBUG:
    //for(unsigned int j = 0; j < this->num_clusters; j++){ 
    //    std::cout << "Cluster #" << j << ": ";
    //    for(unsigned int i = 0; i < clust_idxs[j].size(); i++){
    //        std::cout << " " << clust_idxs[j][i];
    //    }
    //    std::cout << std::endl;
    //}

    for(int j = 0; j < this->num_clusters; j++){
        std::vector<double> curr_data;
        for(auto &idx : clust_idxs[j])
            curr_data.push_back(this->data[idx]);
        this->unique_values[j].sample_given_data(curr_data);
    }

    // std::cout << std::endl; // DEBUG
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal2<Hierarchy, Hypers, Mixture>::eval_density(
        const std::vector<double> grid){
    this->density.first = grid;

    Eigen::VectorXd dens(grid.size());
    double M = this->mixture.get_totalmass();
    int n = this->data.size();
    IterationOutput state;

    for(unsigned int iter = 0; iter < this->chain.state_size(); iter++){
        // for each iteration of the algorithm

        state = *(this->chain.mutable_state(iter)); // TODO does it work?
        std::vector<unsigned int> card(state.phi_size(),
            0); // TODO salviamoci ste card da qualche parte
        std::vector<double> params(state.phi(0).params_size());
        for(unsigned int j = 0; j < n; j++){
            card[ state.allocations(j) ] += 1;
        }
        Hierarchy<Hypers> temp_hier(this->unique_values[0].get_hypers());
        for(unsigned int h = 0; h < state.phi_size(); h++){
            for(int k = 0; k < state.phi(h).params_size(); k++){
                params[k] = state.phi(h).params(k);
            }
            temp_hier.set_state(params);

            dens += card[h] * temp_hier.like(grid) / (M+n);
        }

        // Component from G0 (exploit conjugacy using the explicit expression)
         dens += M * temp_hier.eval_marg(grid) / (M+n); 
    }

    // DEBUG:
    // for(int i = 0; i < grid.size(); i++)
    //     std::cout << dens(i) << " ";
    // std::cout << std::endl;

    this->density.second = dens / this->chain.state_size();

    //DEBUG:
    // for(int i = 0; i < grid.size(); i++)
    //     std::cout << this->density.second(i) << " ";
    // std::cout << std::endl;
}


#endif // NEAL2_IMP_HPP