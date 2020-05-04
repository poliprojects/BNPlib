import numpy as np


class ConditionalMH(object):
    """
    This class implements Aglorithm 7.1 in Moller-Waagpetersen (2003)
    """
    def __init__(self, n_points, papangelou, proposal_rng, proposal_dens):
        self.n_points = n_points
        self.papangelou = papangelou
        self.proposal_rng = proposal_rng
        self.proposal_dens = proposal_dens
        self.state = np.zeros(n_points)

    def run_one(self):
        ind = np.random.choice(self.n_points)
        csi = self.proposal_rng(self.state, ind)

        # compute acceptance ratio
        prop = self.state.copy()
        aux = np.delete(self.state, ind, axis=0)
        prop[ind] = csi
        arate = self.papangelou(csi, aux, log=True) + \
            self.proposal_dens(prop, self.state, ind, log=True) - \
            self.papangelou(self.state[ind, :], aux, log=True) - \
            self.proposal_dens(self.state, self.state, ind, log=True)

        if np.log(np.random.uniform()) < arate:
            self.state = prop

    def run(self, nburn, nsamples, init_state):
        out = np.zeros((nsamples, self.n_points, init_state.shape[1]))
        self.state = init_state

        for i in range(nburn):
            self.run_one()

        for i in range(nsamples):
            self.run_one()
            out[i, :, :] = self.state

        return out



