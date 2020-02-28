#ifndef PITYORMIXTURE_HPP
#define PITYORMIXTURE_HPP

class PitYorMixture {

private:
    double strength;
    double discount;

public:
    ~PitYorMixture() = default;

    PitYorMixture(const double strength,const double discount): strength(strength), discount(discount){
        assert(strength > -discount);
        assert(discount<1 && discount>=0);
    }

    double const get_strength(){return strength;}

    void set_strength(const double strength_){strength = strength_;}

    double const get_discount(){return discount;}

    void set_discount(const double discount_){discount = discount_;}


    double const prob_old(const int card, const unsigned int n ){return (card-discount)/(n-1+strength);}
    
    double const prob_new(const unsigned int n, const unsigned int n_unique){return (strength+discount*n_unique)/(n-1+strength);}
};

#endif // DIRICHLETMIXTURE_HPP
