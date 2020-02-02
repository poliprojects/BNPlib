#ifndef SIMPLEMIXTURE_HPP
#define SIMPLEMIXTURE_HPP

#include "includes_universal.hpp"

class SimpleMixture {

private:
    double totalmass;

public:
    ~SimpleMixture() = default;

    SimpleMixture(double totalmass): totalmass(totalmass) {
      assert(totalmass>=0);
    }

    double const get_totalmass(){return totalmass;}

};

#endif // SIMPLEMIXTURE_HPP
