#include "base_pp.hpp"

BasePP::BasePP(const MatrixXd &ranges): ranges(ranges) {
    dim = ranges.cols();
    diff_range = (ranges.row(1) - ranges.row(0)).transpose();
    vol_range = diff_range.prod();
}

void BasePP::set_ranges(const MatrixXd &ranges) {
    this->ranges = ranges;
    dim = ranges.cols();
    std::cout << "ranges: \n" << ranges << std::endl;
    diff_range = (ranges.row(1) - ranges.row(0)).transpose();
    vol_range = diff_range.prod();

    initialize();
}

MatrixXd BasePP::sample_uniform(int npoints)
{
    MatrixXd out(npoints, dim);
    for (int j=0; j < dim; j++) {
        for (int i=0; i < npoints; i++) {
            out(i, j) = uniform_rng(
                ranges(0, j), ranges(1, j), Rng::Instance().get());
        }
    }    

    return out;
}

void BasePP::sample_given_active(
        const MatrixXd &active, MatrixXd *non_active, double psi_u) {
    int npoints = non_active->rows();
    double c_star_na = c_star * psi_u;

    double rsecond = uniform_rng(0, 1, Rng::Instance().get());
    birth_prob = c_star / (c_star + npoints);
    birth_arate = -1;

    if (rsecond < birth_prob) {
        // BIRTH MOVE
        VectorXd xi = phi_star_rng();
        MatrixXd aux(active.rows() + npoints, dim);
        aux << active, *non_active;
        birth_arate = papangelou(xi, aux) - phi_star_dens(xi) +
          std::log(psi_u);


        double rthird = uniform_rng(0, 1, Rng::Instance().get());
        if (std::log(rthird) < birth_arate)
        {
            non_active->conservativeResize(npoints + 1, dim);
            non_active->row(npoints) = xi;
        }
    } else {
        // Death Move
        if (npoints == 0)
            return;

        VectorXd probas = VectorXd::Ones(npoints) / npoints;
        int ind = categorical_rng(probas, Rng::Instance().get()) - 1;

        delete_row(non_active, ind);

        // if (ind < npoints) {
        //     non_active->block(ind, 0, npoints - ind, dim) = \
        //         non_active->block(ind + 1, 0, npoints - ind, dim);
        // }

        // non_active->conservativeResize(npoints - 1, dim);
    }
    return;
}