#ifndef DIRICHLETMIXTURE_HPP
#define DIRICHLETMIXTURE_HPP

class DirichletMixture {

private:
    double totalmass;

public:
    ~DirichletMixture() = default;

    DirichletMixture(const double &totalmass): totalmass(totalmass){
        assert(totalmass >= 0);
    }

    // Compute probabilities
    double const prob_existing_cluster(const int &card, const unsigned int &n){
    	return card/(n-1+totalmass);
    }
    
    double const prob_new_cluster(const unsigned int &n,
    	const unsigned int &n_unique){
    	return totalmass/(n-1+totalmass);
    }

    // Getters and setters
    double const get_totalmass(){return totalmass;}
    void set_totalmass(const double &totalmass_){totalmass = totalmass_;}
};

#endif // DIRICHLETMIXTURE_HPP
