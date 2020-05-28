#ifndef DIRICHLETMIXTURE_HPP
#define DIRICHLETMIXTURE_HPP


//! Class that represents the Dirichlet process mixture model.

//! This class represents a particular mixture model for iterative BNP algo-
//! rithms, called the Dirichlet process. It represents the distribution of a
//! random probability measure that fulfills a certain property involving the
//! Dirichlet distribution. In terms of the algorithms, it translates to a mix-
//! ture that assigns a weight M, called the total mass parameter, to the crea-
//! tion of a new cluster, and weights of already existing clusters are propor-
//! tional to their cardinalities.

class DirichletMixture {
protected:
    //! Total mass parameters
    double totalmass;

public:
    // DESTRUCTOR AND CONSTRUCTORS
    ~DirichletMixture() = default;
    DirichletMixture() = default;
    DirichletMixture(const double totalmass_): totalmass(totalmass_){
        assert(totalmass >= 0);
    }

    // PROBABILITIES FUNCTIONS (TODO)
    double prob_existing_cluster(const unsigned int card, const unsigned int n)
        const {
        return card/(n-1+totalmass);
    }

    double prob_new_cluster(const unsigned int n, const unsigned int n_clust)
        const {
        return totalmass/(n-1+totalmass);
    }

    // GETTERS AND SETTERS
    double get_totalmass() const {return totalmass;}
    void set_totalmass(const double totalmass_){totalmass = totalmass_;}
};


#endif // DIRICHLETMIXTURE_HPP
