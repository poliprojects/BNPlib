#ifndef SIMPLEMIXTURE_HPP
#define SIMPLEMIXTURE_HPP

class SimpleMixture {

private:
    double totalmass;

public:
    ~SimpleMixture() = default;

    SimpleMixture(const double totalmass): totalmass(totalmass){
        assert(totalmass >= 0);
    }

    double const get_totalmass(){return totalmass;}

    void set_totalmass(const double totalmass_){totalmass = totalmass_;}

};

#endif // SIMPLEMIXTURE_HPP
