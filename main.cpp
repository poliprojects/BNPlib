#include <iostream>
#include <fstream>

//#include <boost/random/random_number_generator.hpp>
//#include <boost/random/detail/qrng_base.hpp>

#include "includes_main.hpp"
#include "math.h"

int main(int argc, char *argv[]){
	
    // 3D-vectorial data
    Eigen::MatrixXd data;
    data << 1.3, 0.9, 8.8, 2.0, -1.3,
            2.3, 5.1, 4.4, 0.0, -3.2,
            3.0, 4,0, 3.3, 1.1, +1.5;

    Eigen::VectorXd mu0(3.0, 3.0, 3.0);
    Eigen::MatrixXd lambda0 = 2 * Matrix<double, 3, 3>::Identity();
    double totalmass = 1.0;
    int n_aux = 3;

    HypersDummy hy(mu0, lambda0);
    DirichletMixture mix(totalmass); // total mass

    Neal8<HierarchyDummy, HypersDummy, DirichletMixture> sampler(
        data, n_aux, mix, hy);
	
    // Run sampler
    sampler.run();

    return 0;
}
