import logging

import numpy as np
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde
from sklearn.metrics import pairwise_distances

import pp_mix.protos.py.params_pb2 as params_pb2
from pp_mix.protos.py.params_pb2 import Params


def make_params(pp_params, prec_params, jump_params):
    params = Params()
    if isinstance(pp_params, params_pb2.StraussParams):
        params.strauss.CopyFrom(pp_params)
    elif isinstance(pp_params, params_pb2.NrepParams):
        params.nrep.CopyFrom(pp_params)
    else:
        logging.error("pp_params type not valid")

    if isinstance(prec_params, params_pb2.WishartParams):
        params.wishart.CopyFrom(prec_params)
    elif isinstance(prec_params, params_pb2.FixedMultiPrecParams):
        params.fixed_multi_prec.CopyFrom(prec_params)
    elif isinstance(prec_params, params_pb2.GammaParams):
        params.gamma_prec.CopyFrom(prec_params)
    elif isinstance(prec_params, params_pb2.FixedUnivPrecParams):
        params.fixed_univ_prec.CopyFrom(prec_params)
    else:
        logging.error("prec_params type not valid")

    if isinstance(jump_params, params_pb2.GammaParams):
        params.gamma_jump.CopyFrom(jump_params)
    else:
        logging.error("jump_params not recognized")

    return params



def make_default_strauss(data, nstar=10, m_max=30):
    params = params_pb2.StraussParams()
    pdist = pairwise_distances(data).reshape((-1, ))
    grid = np.linspace(np.min(pdist), np.max(pdist), 200)
    if len(pdist) > 10000:
        pdist = np.random.choice(pdist, 10000, False)
    dens_estimate = gaussian_kde(pdist).evaluate(grid)

    params.init.R = grid[argrelextrema(dens_estimate, np.less)[0][0]]
    params.init.gamma = np.exp(-nstar)

    ranges = np.vstack([np.min(data, axis=0), np.max(data, axis=0)])
    vol = np.prod(np.diff(ranges, axis=0))
    params.prior.beta_l = 1.0 / (2 * vol)
    params.prior.beta_u = 2 * m_max / vol


    params.init.beta = params.prior.beta_l 

    return params


def check_params(params, data):
    if data.ndim == 1:
        _check_prec_univariate_params(params, data)
    else:
        _check_prec_multivariate_params(params, data)


def _check_prec_univariate_params(params, data):
    if params.WhichOneof("prec_params") not in \
            ("fixed_univ_prec", "gamma_prec"):
        raise ValueError(
            "Found {0} as precision parameter, expected one of: {1}".format(
                params.WhichOneof("prec_params"),
                "[{0}]".format(", ".join(("fixed_univ_prec", "gamma_prec")))
            ))

    if params.WhichOneof("prec_params") == "gamma_prec":
        if params.gamma_prec.alpha <= 0:
            raise ValueError(
                "Parameter gamma_prec.alpha sould be strictly greater than 0, "
                "found gamma_prec.alpha={0} instead".format(
                    params.gamma_prec.alpha))

        if params.gamma_prec.beta <= 0:
            raise ValueError(
                "Parameter gamma_prec.beta sould be strictly greater than 0, "
                "found gamma_prec.beta={0} instead".format(
                    params.gamma_prec.beta))


def _check_prec_multivariate_params(params, data):
    if params.WhichOneof("prec_params") not in \
            ("fixed_multi_prec", "wishart"):
        raise ValueError(
            "Found {0} as precision parameter, expected one of: {1}".format(
                params.WhichOneof("prec_params"),
                "[{0}]".format(", ".join(("fixed_multi_prec", "wishart")))
            ))

    if params.WhichOneof("prec_params") == "wishart":
        if params.wishart.nu < data.ndim + 1:
            raise ValueError(
                """Parameter wishart.nu sould be strictly greater than {0} + 1,
                 found wishart.nu={1} instead""".format(
                     data.ndim, params.wishart.nu))

        if params.wishart.identity is False:
            raise ValueError(
                "Only 'True' is supported for parametr wishart.identity")

        if params.wishart.dim != data.ndim:
            raise ValueError(
                "Parameter wishart.dim should match the dimension of the data, "
                "found wishart.dim={0}, data.ndim={1}".format(
                    params.wishart.dim, data.ndim))

        if params.wishart.HasField("sigma") and params.wishart.sigma <= 0:
            raise ValueError(
                "Parameter wishart.sigma should be grater than 0, "
                "found wishart.wigma={0} instead".format(params.wishart.sigma))
