#ifndef PITYORMIXTURE_HPP
#define PITYORMIXTURE_HPP

class PitYorMixture {
protected:
    double strength;
    double discount;

public:
    // Destructor and constructor
    ~PitYorMixture() = default;
    PitYorMixture() = default;
    PitYorMixture(const double strength_, const double discount_):
        strength(strength_), discount(discount_){
        assert(strength > -discount);
        assert(0 <= discount && discount < 1);
    }

    // Compute probabilities
    double prob_existing_cluster(const unsigned int card, const unsigned int n)
        const {
        return (card-discount)/(n-1+strength);
    }

    double prob_new_cluster(const unsigned int n,
        const unsigned int n_unique) const {
        return (strength+discount*n_unique)/(n-1+strength);
    }
    
    // Getters and setters
    double get_strength() const {return strength;}
    double get_discount() const {return discount;}
    void set_strength(const double strength_){strength = strength_;}
    void set_discount(const double discount_){discount = discount_;}
};


#endif // DIRICHLETMIXTURE_HPP
