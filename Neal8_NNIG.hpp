#include <tuple>
#include <stan/math/prim/mat.hpp>
#include <type_traits>

#include "includes.hpp"
#include <math.h>   
// N-NIG model == gaussian kernel + N-IG base measure:
// f ~ N(mu,sig^2)
// (mu,sig^2) ~ G
// G ~ DP(M, G0)  with G0 = N-IG

template<template <class> class Hierarchy, class Hypers, class Mixture>
class Neal8{
private:
    unsigned int n_aux=3;
    unsigned int maxiter = 10000; // TODO LATER
    unsigned int burnin = 0;
    std::mt19937 rng;
    int numClusters;
    Mixture mixture;
    //Hypers hy;

    std::vector<data_t> data;
    std::vector<unsigned int> allocations; // the c vector
    std::vector<Hierarchy<Hypers>> unique_values;
    std::vector<Hierarchy<Hypers>> aux_unique_values;



    void initalize(){
        std::default_random_engine generator;
        std::uniform_int_distribution<int> distribution(0,numClusters);
    
        for (int h = 0; h < numClusters; h++) {
          allocations.push_back(h);
        }
   
        for (int j = numClusters; j < data.size(); j++) {
            int num = distribution(generator); //TODO da stan?
            allocations[j] = num;
        }
    }

    void step(){
        sample_allocations();
        sample_unique_values();
    }

    void sample_allocations(){
        // Other ideas:
        // * our own for loop for k and bool (ci is a singleton)
        // * function from std count distinct values in vector
        // * using a (multi)map?

        // Initialize some relevant variables
        unsigned int k, n_unique, singleton;
        unsigned int n=data.size();
    
        

        
 
        for(int i=0; i<n; i++){ // for each data unit data[i]

            // Initialize cardinalities of unique values
            std::vector<int> card(unique_values.size(), 0);      
            for(int j=0; j<n; j++)
                card[ allocations[j] ] += 1;
    
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
            for(int j=singleton; j<n_aux; j++){
                aux_unique_values[j].draw();
            }

            // Draw a NEW value for ci
            Eigen::MatrixXd probas(n_unique+n_aux,1); //k or n_unique
            //Matrix<double, Dynamic, 1> VectorXd

            auto M = mixture.get_totalmass();
            double tot=0.0;

            for(int k=0; k<n_unique ; k++){ // if datum i is a singleton, then
                // card[k] when k=allocations[i] is equal to 0 -> probas[k]=0
    
                // TODO LATER "meglio in logscale" (?)
                probas(k,0) = card[k] * unique_values[k].log_like(data[i]) / (
                    n-1+M);

                tot+=probas(k,0);
            }

            for(int k=0; k<n_aux; k++){
                probas(n_unique+k,0) = (M/n_aux) *
                    aux_unique_values[k].log_like(data[i]) / (n-1+M);
                tot += probas(n_unique+k,0);
               }
            probas = probas * (1/tot);
      
            for(int i=0; i<probas.size(); i++){
                std::cout << "probas_" << probas(i,0) << std::endl; // DEBUG
            }
 


            unsigned int c_new = stan::math::categorical_rng(probas, rng) -1;
            
            std::cout<<"c_new: "<<c_new<<std::endl; // DEBUG

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
                        if (c > tmp)
                            c -= 1;
                    }
                } // end of else

            } // end of if(singleton == 1)

            else{ // if singleton == 0

                if (c_new>=n_unique){ // case 3 of 4: NOT SINGLETON - AUX
                    unique_values.push_back(aux_unique_values[c_new-n_unique]);
                    card.push_back(1);
                    allocations[i] = n_unique;

                }
                else{ // case 4 of 4: NOT SINGLETON - OLD VALUES
                    allocations[i] = c_new;
                    card[c_new] += 1;
                }
            } // end of else

        } // end of for(int i=0; i<n; i++) loop

    } // end of sample_allocations()



    void sample_unique_values(){

        numClusters=unique_values.size();

        std::vector<std::vector<unsigned int>> clust_idxs(numClusters);
        unsigned int n = allocations.size();


        for(unsigned int i=0; i<n; i++){ // save different cluster in each row
            clust_idxs[ allocations[i] ].push_back(i);
        }


        // DEBUG:
        for(int j=0; j<numClusters; j++){ 
            std::cout << "Cluster #" << j << ": ";

            for (unsigned int i=0; i<clust_idxs[j].size(); i++)
                std::cout << " " << clust_idxs[j][i];

            std::cout << std::endl;
        }

        for (unsigned int j=0; j< numClusters; j++) {
            std::vector<data_t> curr_data;

            for ( auto &idx : clust_idxs[j] )
                curr_data.push_back( data[idx] );

            unique_values[j].sample_given_data(curr_data);
        }


    std::cout << std::endl; // DEBUG

    }



    void save_iteration(unsigned int iter){ // TODO LATER
        std::cout << "Iteration n. " << iter << " / " << maxiter << std::endl;
        print();
    }

    void print() {
        for (int h = 0; h < numClusters; h++) {
            std::cout << "Cluster # " << h << std::endl;
            std::cout << "Parameters: ";

            for (auto c : unique_values[h].get_state()){
                std::cout << c << " " << std::endl;
            }
            std::cout << std::endl;
        }
    }


public:
    // Constructors and destructors:
    ~Neal8() = default;
    Neal8(const std::vector<data_t> & data, int numClusters,int n_aux,
        const Mixture & mix,const Hypers &hy):
        data(data), numClusters(numClusters), n_aux(n_aux), mixture(mix) {
            Hierarchy<Hypers> hierarchy(std::make_shared<Hypers> (hy));
            for (int h=0; h<numClusters; h++) {
                unique_values.push_back(hierarchy);
            }
            for (int h=0; h<n_aux; h++) {
                aux_unique_values.push_back(hierarchy);
            }
    }
    // If no # initial clusters is given, it will be equal to the data size:
    Neal8(std::vector<data_t> &data, int n_aux, const Mixture & mix,
        const Hypers &hy): Neal8(data, data.size(), n_aux, mix, hy) {}

    // Running tool
    void run(){
        initalize();
        unsigned int iter = 0;
        while(iter < maxiter){
            step();    
            if(iter >= burnin)
              save_iteration(iter);
            iter++;
        }
    }

}; // end of Class Neal8
