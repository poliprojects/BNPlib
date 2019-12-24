// N-NIG model == gaussian kernel + N-IG base measure:
// f ~ N(mu,sig^2)
// (mu,sig^2) ~ G
// G ~ DP(M, G0)  with G0 ~ N-IG


// Normal likelihoood, Normal Inverse Gamma hierarchy
template<class Hypers> //Hypers = TupleWrapper, distro, ...
class NNIGHierarchy {
protected:
    using std::tuple<par_t, par_t, par_t, par_t> = state_tuple_t;

    unsigned int rng = 20191225;
    state_tuple_t state; // current values for G0's parameters: in order, mu_0,
                         // sig2_0, alpha, beta

    std::shared_ptr<Hypers> hypers;

public:
    // Contructors:
    // TODO!

    // Getters/setters:
    state_tuple_t get_state(){return state;}
    void set_state(const state_tuple_t &s){state = s;}
    void set_state(int pos, par_t val){state(pos) = val;}

    double log_like(data_t datum) {
        return stan::math::normal_lpdf(datum, state(0), state(1));
    }

    void draw() {
        real sigmaNew = stan::math::inv_gamma_rng(hypers.get_alpha0(),
            hypers.get_beta0(), rng);
        real muNew = stan::math::normal_rng(hypers.get_m0(), hypers.get_sig20(),
            rng);
        state(0) = muNew;
        state(1) = sigmaNew;
        }

    void sample_given_data(std::vector<data_t> data) { // TODO
        // Get current values of parameters
        auto mu0    = hypers.get_mu0();
        auto sig20  = hypers.get_sig20();
        auto alpha0 = hypers.get_alpha0();
        auto beta0  = hypers.get_beta0();

        // Compute posterior parameters // TODO
        auto mu_post = ...;
        auto sig_post = ...;
        auto alpha_post = ...;
        auto beta_post = ...;

        // Get a sample
        par_t sigma_new = stan::math::inv_gamma_rng(alpha_post, beta_post, rng);
        par_t mu_new = stan::math::normal_rng(mu_post, sig_post*sigma_new, rng);
        state(0) = mu_new;
        state(1) = sigma_new;
    }

}




template<class Hierarchy, class Hypers, unsigned int n_aux = 3>
class Neal8{
private:
    unsigned int m = n_aux; // TODO
    unsigned int maxiter = 1000;
    unsigned int burnin = 0;

    //std::vector<Hierarchy> hierarchies;
    std::vector<data_t> data;
    std::array<n, unsigned int> allocations; // the c vector
    std::vector<Hierarchy> unique_values;
    std::array<m, Hierarchy> aux_unique_values;
    
    Hypers hypers;


    void initalize(){
        for(int i=0; i<n; i++)
            allocations[i] = i; // one datum per cluster
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

        // Initialize cardinalities of unique values
        std::vector<unsigned int> card(unique_values.size(), 0);
        for(int j=0; j<n; j++)
            card[ allocations[j] ] += 1;

        for(int i=0; i<n; i++){ // for each data unit data[i]

            singleton = 0;
            n_unique = unique_values.size();

            if(card[ allocations[i] ] == 1){ // datum i is a singleton
                k = n_unique - 1;
                aux_unique_values[0].set_state( unique_values[allocations[i]
                    ].get_state() ); // move phi value in aux
                singleton = 1;
            }
            else{
                k = n_unique;
            }

            card[ allocations[i] ] -= 1;


            // Draw the aux from G0
            for(int j=singleton; j<m; j++){
                aux_unique_values[j].draw();
            }

            // Draw a NEW value for ci
            Eigen::Matrix<double, k+m, 1> probas;

            for(int k=0; k<n_unique ; k++){ // if datum i is a singleton, then
                // card[k] when k=allocations[i] is equal to 0 -> probas[k]=0

                auto M = mixture.get_totalmass();

                // TODO giusto? "meglio in logscale" (?)
                probas(k) = card[k] * unique_values[k].loglike(data[i]) / (
                    n-1+M);
            }
            for(int k=0; k<m ; k++){
                probas(n_unique+k) = (M/m) *
                    aux_unique_values[k].loglike(data[i]) / (n-1+M);
            }

            unsigned int c_new = stan::math::categorical_rng(probas,rng);
            if(singleton == 1){
                if (c_new >= n_unique){ // case 1 of 4: SINGLETON - AUX
                    unique_values[ allocations[i] ].set_state(
                        aux_unique_values[c_new-n_unique].get_state());
                    card[ allocations[i] ] += 1;
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
                    unique_values.push_back(aux_unique_values[c_new-n_unique]);
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


    void save_iteration(unsigned int iter){ // TODO
        std::cout << "Iteration n. " << iter << " / " << maxiter << std::endl;
    }


public:
    // Constructors
    // TODO!

    void run(){
        initalize();
        unsigned int iter = 0;
        while(iter < maxiter){
            step();
            if(iter >= burnin)
                save_iteration(iter);
        }
    }

} // end of Class Neal8



int main() {
    Neal8<NNIGHierarchy<HypersFixed>> sampler;
}
