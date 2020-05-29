#ifndef BASEMIXTURE_HPP
#define BASEMIXTURE_HPP


//! Abstract base class for a generic mixture model

//! This class represents a mixture model object to be used in a BNP iterative
//! algorithm. By definition, a mixture is a probability distribution that inte-
//! grates over a density kernel to generate the actual distribution for the da-
//! ta. However, in the context of this library, where a clustering structure is
//! generated on the data, a certain mixture translates to a certain way of
//! weighing the insertion of data in old clusters vs the creation of new clus-
//! ters. Therefore any mixture object inheriting from the class must have me-
//! thods that provide the probabilities for the two aforementioned events. The
//! class will then have its own parameters, and maybe even prior distributions
//! on them.

class BaseMixture {
public:
    // DESTRUCTOR AND CONSTRUCTORS
    virtual ~BaseMixture() = default;
    BaseMixture() = default;


    // PROBABILITIES FUNCTIONS
    //! Mass probability for choosing an already existing cluster

    //! \param card Cardinality of the cluster
    //! \param n    Total number of data points
    //! \return     Probability value
    virtual double mass_existing_cluster(const unsigned int card,
        const unsigned int n) const = 0;


    //! Mass probability for choosing a newly created cluster

    //! \param n_clust Number of clusters
    //! \param n       Total number of data points
    //! \return        Probability value
    virtual double mass_new_cluster(const unsigned int n_clust,
        const unsigned int n) const = 0;
};


#endif // BASEMIXTURE_HPP
