#ifndef PITYORMIXTURE_HPP
#define PITYORMIXTURE_HPP


//! Class that represents the Pitman-Yor process mixture model.

//! This class represents a particular mixture model for iterative BNP algo-
//! rithms, called the Pitman-Yor process. It has two parameters, the strength
//! and the discount. It is a generalized version of the Dirichlet process,
//! which has discount = 0 and strength = total mass. In terms of the algo-
//! rithms, it translates to a mixture that assigns a weight to the creation of
//! a new cluster proportional to their cardinalities, but reduced by the dis-
//! count factor, while the weight for a newly created cluster is the remaining
//! one counting the total amount as the sample size increased by the strength.

class PitYorMixture {
protected:
    //! Strength and discount parameters
    double strength, discount;

public:
    // DESTRUCTOR AND CONSTRUCTORS
    ~PitYorMixture() = default;
    PitYorMixture() = default;
    PitYorMixture(const double strength_, const double discount_):
        strength(strength_), discount(discount_){
        assert(strength > -discount);
        assert(0 <= discount && discount < 1);
    }

    // PROBABILITIES FUNCTIONS (TODO)
    double mass_existing_cluster(const unsigned int card, const unsigned int n)
        const {
        return (card-discount)/(n+strength);
    }

    double mass_new_cluster(const unsigned int n_clust, const unsigned int n)
    	const {
        return (strength+discount*n_clust)/(n+strength);
    }

    // GETTERS AND SETTERS
    double get_strength() const {return strength;}
    double get_discount() const {return discount;}
    void set_strength(const double strength_){strength = strength_;}
    void set_discount(const double discount_){discount = discount_;}
};


#endif // PITYORMIXTURE_HPP
