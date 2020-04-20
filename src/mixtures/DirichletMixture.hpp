#ifndef DIRICHLETMIXTURE_HPP
#define DIRICHLETMIXTURE_HPP

class DirichletMixture {
private:
    double totalmass;

public:
    // Destructor and constructor
    ~DirichletMixture() = default;

    DirichletMixture(const double totalmass): totalmass(totalmass){
        assert(totalmass >= 0);
    }

    // Compute probabilities
    double prob_existing_cluster(const unsigned int card, const unsigned int n)
        const {
    	return card/(n-1+totalmass);
    }
    
    double prob_new_cluster(const unsigned int n,
    	const unsigned int n_unique) const {
    	return totalmass/(n-1+totalmass);
    }

    // Getters and setters
    double get_totalmass() const {return totalmass;}
    void set_totalmass(const double totalmass_){totalmass = totalmass_;}
};


#endif // DIRICHLETMIXTURE_HPP
