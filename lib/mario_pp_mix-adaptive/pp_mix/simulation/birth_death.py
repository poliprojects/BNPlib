import jax.numpy as np
from jax import random

import pp_mix.utils as utils


class BirthDeathMH(object):
    """
    This cass implements Aglorithm 7.5
    """
    def __init__(
            self, papangelou, birth_proposal_rng, birth_proposal_dens,
            update_proposal_rng, update_proposal_dens):
        self.papangelou = papangelou
        self.birth_proposal_rng = birth_proposal_rng
        self.birth_proposal_dens = birth_proposal_dens
        self.update_proposal_rng = update_proposal_rng
        self.update_proposal_dens = update_proposal_dens
        self.pbirth = 0.9
        self.state = np.zeros((10, 2))

    def birth_move(self, rng):
        csi = self.birth_proposal_rng(self.state)
        prop = np.vstack([self.state, [csi]])
        arate = self.papangelou(csi, self.state) + (1 - self.pbirth) - \
            np.log(prop.shape[0]) - \
            self.pbirth - self.birth_proposal_dens(csi, self.state)
        
        if np.log(random.uniform(rng.get())) < arate:
            self.state = prop

    def death_move(self, rng):
        npoints = self.state.shape[0]
        
        ind = random.randint(rng.get(), (1,), 0, npoints)
        prop = utils.delete(self.state.copy(), ind, axis=0)

        arate = self.pbirth + \
            self.birth_proposal_dens(self.state[ind, :], prop) - \
            self.papangelou(self.state[ind, :], self.state) - \
            (1 - self.pbirth) + np.log(npoints)

        if np.log(random.uniform(rng.get())) < arate:
            self.state = prop.copy()

    def update_move(self, rng):
        npoints = self.state.shape[0]
        ind = random.randint(rng.get(), (1,), 0, npoints)
        csi = self.update_proposal_rng(self.state, ind)

        # compute acceptance ratio
        prop = self.state.copy()
        aux = utils.delete(self.state, ind, axis=0)
        prop[ind] = csi
        arate = self.papangelou(csi, aux, log=True) + \
            self.update_proposal_dens(prop, self.state, ind, log=True) - \
            self.papangelou(self.state[ind, :], aux, log=True) - \
            self.update_proposal_dens(self.state, self.state, ind, log=True)

        if np.log(random.uniform(rng.get())) < arate:
            self.state = prop

    def run_one(self, q, rng):
        if random.uniform(rng.get()) < q and self.state.shape[0] > 1:
            self.update_move(rng)
        else:
            if random.uniform(rng.get()) < self.pbirth:
                self.birth_move(rng)
            elif self.state.shape[0] > 1:
                self.death_move(rng)

    def run(self, rng, nburn, nsamples, init_state, q=0.5):
        out = [None] * nsamples
        self.state = init_state

        for i in range(nburn):
            self.run_one(q, rng)

        for i in range(nsamples):
            self.run_one(q, rng)
            out[i] = self.state

        return out


class SpatialBirthAndDeath(object):
    """
    This class implements Algorithm 11.3
    """
    def __init__(self, dens, papangelou, phi_star_rng, phi_star_dens, c_star):
        self.jump_times = None
        self.y = None # realizations
        self.x = None # realizations from the dominating process
        self.T = None
        self.papangelou = papangelou
        self.phi_star_rng = phi_star_rng
        self.phi_star_dens = phi_star_dens
        self.c_star = c_star

    def initialize(self, init_state):
        self.y = init_state
        self.x = init_state
        self.T = 0
        self.npoints = init_state.shape[0]

    def run(self, rng, nburn, nsamples, init_state):
        out = [None] * nsamples
        self.initialize(init_state)

        for i in range(nburn):
            self.run_one(rng)

        for i in range(nsamples):
            self.run_one(rng)
            out[i] = self.y

        return out

    def run_one(self, rng):
        rpime = random.uniform(rng.get())
        rsecond = random.uniform(rng.get())
        # self.T += np.log(-rpime) / (self.c_star + self.npoints)

        if rsecond <= self.c_star / (self.c_star + self.npoints):
            # BIRTH MOVE
            csi = self.phi_star_rng()
            self.x = np.vstack([self.x, csi])

            arate = self.papangelou(csi, self.y) - self.phi_star_dens(csi)
            if np.log(random.uniform(rng.get())) < arate:
                self.y = np.vstack([self.y, csi])
                self.npoints += 1
        else:
            # DEATH MOVE
            if self.npoints == 0:
                return

            ind = random.randint(rng.get(), (1,), 0, self.npoints)
            deleted = self.x[ind, :]
            self.x = utils.delete(self.x, ind, axis=0)
            wh = np.where(self.y == deleted)[0]
            if len(wh):
                self.y = utils.delete(self.y, wh, axis=0)
                self.npoints -= 1
            


