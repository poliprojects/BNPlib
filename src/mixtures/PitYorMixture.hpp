#ifndef PITYORMIXTURE_HPP
#define PITYORMIXTURE_HPP

class PitYorMixture {

private:
    double strength;
    double discount;

public:
    ~PitYorMixture() = default;

    PitYorMixture(const double strength, const double discount):
        strength(strength), discount(discount){
        assert(strength > -discount);
        assert(0 <= discount && discount<1);
    }

    // Compute probabilities
    double const prob_existing_cluster(const int card, const unsigned int n){
        return (card-discount)/(n-1+strength);
    }

    double const prob_new_cluster(const unsigned int n,
        const unsigned int n_unique){
        return (strength+discount*n_unique)/(n-1+strength);
    }
    
    // Getters and setters
    double const get_strength(){return strength;}
    double const get_discount(){return discount;}
    void set_strength(const double strength_){strength = strength_;}
    void set_discount(const double discount_){discount = discount_;}
};

#endif // DIRICHLETMIXTURE_HPP
