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
