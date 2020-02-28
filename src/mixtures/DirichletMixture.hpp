#ifndef DIRICHLETMIXTURE_HPP
#define DIRICHLETMIXTURE_HPP

class DirichletMixture {

private:
    double totalmass;

public:
    ~DirichletMixture() = default;

    DirichletMixture(const double totalmass): totalmass(totalmass){
        assert(totalmass >= 0);
    }

    double const get_totalmass(){return totalmass;}

    void set_totalmass(const double totalmass_){totalmass = totalmass_;}

    double const prob_old(const int card, const unsigned int n ){return card/(n-1+totalmass);}
    
    double const prob_new(const unsigned int n, const unsigned int n_unique){return totalmass/(n-1+totalmass);}
};

#endif // DIRICHLETMIXTURE_HPP
