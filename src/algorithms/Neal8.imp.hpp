#ifndef NEAL8_IMP_HPP
#define NEAL8_IMP_HPP

#include "Neal8.hpp"

template<template <class> class Hierarchy, class Hypers, class Mixture>
const void Neal8<Hierarchy, Hypers, Mixture>::print_startup_message(){
    std::cout << "Running Neal8 algorithm (with m=" << n_aux <<
        " auxiliary blocks)..." << std::endl;
}


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal8<Hierarchy, Hypers, Mixture>::sample_allocations(){
    // TODO Other ideas:
    // * our own for loop for k and bool (ci is a singleton)
    // * function from std count distinct values in vector
    // * using a (multi)map?
    // Initialize some relevant variables
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
            k = n_unique - 1;
            aux_unique_values[0].set_state( unique_values[ allocations[i]
                ].get_state() ); // move phi value in aux
            singleton = 1;
        }
        else{
            k = n_unique;
        }

        card[ allocations[i] ] -= 1;

        // Draw the aux from G0
        for(int j = singleton; j < n_aux; j++){
            aux_unique_values[j].draw();
        }

        // Draw a NEW value for ci
        Eigen::VectorXd probas(n_unique+n_aux); //k or n_unique

        //auto M = mixture.get_totalmass();
        double tot = 0.0;
        for(int k = 0; k < n_unique ; k++){ // if datum i is a singleton, then
            // card[k] when k=allocations[i] is equal to 0 -> probas[k]=0
            probas(k) = mixture.prob_existing_cluster(card[k],n) *
            	unique_values[k].like(data[i]);
            tot += probas(k);
        }

        for(int k = 0; k < n_aux; k++){
            probas(n_unique+k) = mixture.prob_new_cluster(n, n_unique) *
                aux_unique_values[k].like(data[i]) / n_aux;
            tot += probas(n_unique+k,0);
        }
        probas = probas / tot;

        //for(int i = 0; i < probas.size(); i++){
        //    std::cout << "probas_" << probas(i) << std::endl; // DEBUG
        //}
        unsigned int c_new = stan::math::categorical_rng(probas, rng) - 1;

        //std::cout<<"c_new: "<<c_new<<std::endl; // DEBUG

        if(singleton == 1){
            if(c_new >= n_unique){ // case 1 of 4: SINGLETON - AUX
                unique_values[ allocations[i] ].set_state(
                    aux_unique_values[c_new-n_unique].get_state());
                card[ allocations[i] ] += 1;
            }
            else{ // case 2 of 4: SINGLETON - OLD VALUE
                unique_values.erase(
                    unique_values.begin()+allocations[i] );
                card.erase( card.begin()+allocations[i] );
                card[c_new] += 1;
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
            if(c_new >= n_unique){ // case 3 of 4: NOT SINGLETON - AUX
                unique_values.push_back(aux_unique_values[c_new-n_unique]);
                card.push_back(1);
                allocations[i] = n_unique;
            }
            else{ // case 4 of 4: NOT SINGLETON - OLD VALUES
                allocations[i] = c_new;
                card[c_new] += 1;
            }
        } // end of else

    } // end of for(int i = 0; i < n; i++) loop

} // end of sample_allocations()


template<template <class> class Hierarchy, class Hypers, class Mixture>
void Neal8<Hierarchy, Hypers, Mixture>::eval_density(
        const std::vector<double> grid){
    density.first = grid;

    //std::ofstream file;
    //unsigned int step = 1000;
    //file.open("dens_estimate_iterations.csv");

    Eigen::VectorXd dens(grid.size());
    double M = mixture.get_totalmass();
    int n = data.size();
    IterationOutput state;

    for(int iter = 0; iter < chain.state_size(); iter++){
        // for each iteration of the algorithm
        //std::cout << iter << std::endl; // DEBUG

        state = *chain.mutable_state(iter);
        std::vector<unsigned int> card(state.phi_size(),
            0); // TODO salviamoci ste card da qualche parte
        std::vector<double> params(state.phi(0).params_size());
        Eigen::VectorXd dens_addendum = Eigen::MatrixXd::Zero(grid.size(), 1);

        for(int j = 0; j < n; j++){
            card[ state.allocations(j) ] += 1;
        }
        Hierarchy<Hypers> temp_hier(unique_values[0].get_hypers());
        for(int h = 0; h < state.phi_size(); h++){
            for(int k = 0; k < state.phi(h).params_size(); k++){
                params[k] = state.phi(h).params(k);
            }
            temp_hier.set_state(params);

            dens_addendum += card[h] * temp_hier.like(grid) / (M+n);
            
        }
    
        // Component from G0
        for(int h = 0; h < n_aux; h++){
            temp_hier.draw();
            dens_addendum += (M/n_aux) * temp_hier.like(grid) / (M+n);
        }

    dens += dens_addendum;

        //if(iter % step == 0){
           // for(int i=0; i<dens_addendum.size()-1; i++){

         //   file << dens_addendum(i)<< ",";
        //    }
        //    file <<dens_addendum(dens_addendum.size()-1) << std::endl;
        //}
    }

    // DEBUG:
    // for(int i = 0; i < grid.size(); i++)
    //     std::cout << dens(i) << " ";
    // std::cout << std::endl;

    density.second = dens / chain.state_size();

    //DEBUG:
    // for(int i = 0; i < grid.size(); i++)
    //     std::cout << density.second(i) << " ";
    // std::cout << std::endl;

    //file.close();
}


#endif // NEAL8_IMP_HPP
