import numpy as np
from scipy.linalg import cho_factor


class PrecMat(object):
    def __init__(self, prec):
        self.prec = prec
        self.prec_cho = cho_factor(prec)[0]
        self.log_det_inv = - 2 * np.sum(np.log(self.prec_cho.diagonal()))
