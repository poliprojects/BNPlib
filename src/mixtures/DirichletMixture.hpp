#ifndef DIRICHLETMIXTURE_HPP
#define DIRICHLETMIXTURE_HPP

#include "BaseMixture.hpp"


//! Class that represents the Dirichlet process mixture model.

//! This class represents a particular mixture model for iterative BNP algo-
//! rithms, called the Dirichlet process. It represents the distribution of a
//! random probability measure that fulfills a certain property involving the
//! Dirichlet distribution. In terms of the algorithms, it translates to a mix-
//! ture that assigns a weight M, called the total mass parameter, to the crea-
//! tion of a new cluster, and weights of already existing clusters are propor-
//! tional to their cardinalities.

class DirichletMixture : public BaseMixture {
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


    // PROBABILITIES FUNCTIONS
    //! Mass probability for choosing an already existing cluster

    //! \param card Cardinality of the cluster
    //! \param n    Total number of data points
    //! \return     Probability value
    double mass_existing_cluster(const unsigned int card, const unsigned int n)
        const override {
        return card/(n+totalmass);
    }


    //! Mass probability for choosing a newly created cluster

    //! \param n_clust Number of clusters
    //! \param n       Total number of data points
    //! \return        Probability value
    double mass_new_cluster(const unsigned int n_clust, const unsigned int n)
        const override {
        return totalmass/(n+totalmass);
    }

    // GETTERS AND SETTERS
    double get_totalmass() const {return totalmass;}
    void set_totalmass(const double totalmass_){totalmass = totalmass_;}
};


#endif // DIRICHLETMIXTURE_HPP
