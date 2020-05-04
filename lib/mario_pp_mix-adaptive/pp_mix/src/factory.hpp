#ifndef FACTORY_HPP
#define FACTORY_HPP

#include <memory>

#include "precs/base_prec.hpp"
#include "precs/fixed_prec.hpp"
#include "precs/wishart.hpp"
#include "precs/gamma.hpp"

#include "jumps/base_jump.hpp"
#include "jumps/gamma.hpp"

#include "point_process/base_pp.hpp"
#include "point_process/strauss_pp.hpp"
#include "point_process/nrep_pp.hpp"

#include "conditional_mcmc.hpp"

#include "../protos/cpp/params.pb.h"

BasePP* make_pp(const Params& params);

BasePP* make_strauss(const StraussParams& params);

BasePP* make_nrep(const NrepParams&params);


BaseJump* make_jump(const Params& params);

BaseJump* make_gamma_jump(const GammaParams& params);


BasePrec* make_prec(const Params& params);

BasePrec *make_fixed_prec(const FixedMultiPrecParams &params);

BasePrec* make_wishart(const WishartParams& params);

BasePrec* make_fixed_prec(const FixedUnivPrecParams& params);

BasePrec* make_gamma_prec(const GammaParams& params);

#endif