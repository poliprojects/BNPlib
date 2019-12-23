// Modello N-NIG: kernel gaussiano, G0 N-IG:
// f ~ N(mu,sig^2)
// (mu,sig^2) ~ G
// G ~ DP(M, G0)  with G0 ~ N-IG


// Normal likelihoood, Normal Inverse Gamma Hierarchy
template<typename Hypers> //Hypers = TupleWrapper, distro, ...
class NNIGHierarchy {
protected:
    std::tuple<real, real> state;
    //std::tuple<> params; // mu_0, k, alpha, beta


    //qualcosa che mi rappresenta i parametri ma possono cambiare
    Hypers* hypers;
public:
    // getters, such as:
    std::tuple getState();

    double log_like(real datum) {
        return stan::math::normal_lpdf(datum, state-1, state-2);
    }

    void draw() {
        real sigmaNew = stan::math::inv_gamma_rng(hypers.get_alpha0(),hypers.get_beta0(), rng);
        real muNew = stan::math::normal_rng(hypers.get_m0(), hypers.get_k0(), rng);
        state(0) = muNew;
        state(1) = sigmaNew;
        }


    void sample(std::vector<real> data) {
        // prendo mu_0, k0, alpha0, beta0 da hypers:
        mu_0 = hypers.current_value // ?

        ... // calcolo mu_post, k_post, alpha_post, beta_post

        real sigmaNew = stan::math::inv_gamma_rng(alpha_post, beta_post, rng);
        real muNew = stan::math::normal_rng(mu_post, kpost * sigmaNew, rng);
        state(0) = muNew;
        state(1) = sigmaNew;
        }



}



template<Hierarchy, Hypers>
class Neal8{
private:
    int m=3; // TODO
    // std::vector<Hierarchy> hierarchies;
    std::vector<real> data;
    std::array<n,int> allocations; // the c vector
    //std::vector<param_t> unique_values;
    //std::array<m,param_t> aux_unique_values;

    std::vector<Hierarchy> unique_values;
    std::array<m,Hierarchy> aux_unique_values;

    Hypers hypers; // see Hypers class in other files

    // Idea by Mario
    void sample_c() {
        for (i=0; i<numData; i++) {
            old = c[i];
            probas = Eigen::VectorXd(numClusters);
            for (c = 0; c < numClusters; c++) {

                // meglio in logscale
                probas(c) = dataPerClusters[c] * hierarchies[c].like(data[i]);
            }

            probas;

            c[i] = stan::math::categorical_rng(probas, rng);
            dataPerClusters[old] -= 1;
            dataPerClusters[c[i]] += 1;
        }
    }





    void run(){
        initalize();
        while(!converged){
            step();
        }
    }

    void initalize(){
        //...
        return;
    }

    void step(){
        sample_allocations();
        sample_unique_values();
    }


    void sample_allocations(){
        // Other ideas:
        // * ideas: our own for loop for k and bool (ci is a singleton)
        // * ideas: function from std count distinct values in vector
        // * ideas: using a (multi)map?

        // Initialize some relevant variables
        unsigned int k, n_unique, singleton;

        // Initialize cardinalities of unique values
        std::vector<unsigned int> card(unique_values.size(), 0);
        for(int j=0; j<n; j++)
            card[ allocations[j] ] += 1;

        for(int i=0; i<n; i++){ // for each data unit data[i]

            singleton = 0;
            n_unique = unique_values.size();

            if(card[ allocations[i] ] == 1){ // datum i is a singleton
                k = n_unique - 1;
                aux_unique_values[0].setState(unique_values[allocations[i]].getState()); // move phi value in aux
                singleton = 1;
            }
            else{
                k = n_unique;
            }

            card[ allocations[i] ] -= 1;


             // draw from G0 the aux

            for(int j=singleton;j<m; j++){
                aux_unique_values[j].draw();
            }


            //draw a NEW value for ci

            auto probas = Eigen::VectorXd(k+m); // or vec initialized to 0

            for(int k=0; k<n_unique ; k++){ // if datum i is a singleton, then
                // card[k] when k=allocations[i] is equal to 0 ->probas[k]=0
                probas(k) = card[k] * unique_values[k].loglike(data[i]) / (
                    n-1+mixture.get_totalmass() );
            }
            for(int k=0; k<m ; k++){
                probas(n_unique+k) = (mixture.get_totalmass()/m) *
                    aux_unique_values[k].loglike(data[i])/ ( n-1+mixture.get_totalmass() );
            }

            unsigned int c_new = stan::math::categorical_rng(probas,rng);
            if(singleton == 1){
                if (c_new >= n_unique){ // case 1 of 4: SINGLETON - AUX
                    unique_values[ allocations[i] ].setState(aux_unique_values[c_new-n_unique].getState());
                    card[allocations[i]] += 1;
                }
                else{ // case 2 of 4: SINGLETON - OLD VALUE
                    unique_values.erase(
                        unique_values.begin()+allocations[i]-1 );

                    card.erase( card.begin()+allocations[i]-1 );
                    card[c_new] += 1;

                    int tmp = allocations[i];
                    allocations[i] = c_new;
                    for(auto &c : allocations){ // relabeling
                        if (c > tmp)
                            c -= 1;
                    }
                }
            } // end of if(singleton == 1)

            else{ // if singleton == 0

                if (c_new>=n_unique){ // case 3 of 4: NOT SINGLETON - AUX
                    unique_values.push_back( aux_unique_values[c_new-n_unique] );
                    card.push_back(1);
                    allocations[i] = n_unique + 1;
                }
                else{ // case 4 of 4: NOT SINGLETON - OLD VALUES
                    allocations[i] = c_new;
                    card[c_new] += 1;
                }

            }
        } // end of for(int i=0; i<n; i++) loop

    } // end of sample_allocations()




    void sample_unique_values(){
        std::vector<std::vector<unsigned int>> clust_idxs;
        for(int i=0; i<n; i++) // save different cluster in each row
            clust_idxs[ allocations[i] ].push_back(i);

        for (int j=0; j< clust_idxs.size(); j++) {
            std::vector<data_t> curr_data;

            for ( auto &idx : clust_idxs[j] )
                curr_data.push_back( data[idx] );

            unique_values[j].sample(curr_data);
        }
    }


} // end of Class Neal8



int main() {
    Neal8<NNIGHierarchy<HypersFixed>> sampler;
}
