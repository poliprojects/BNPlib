import json
import numpy as np
import pandas as pd
from collections import defaultdict


def indicatorChainToClusters(rhoChain):
    elems = set(rhoChain.flatten().tolist())
    out = []
    for i in range(rhoChain.shape[0]):
        row = rhoChain[i, ]
        clusters = []
        for e in elems.intersection(row):
            clus = np.where(row == e)[0] + 1
            clusters.append(tuple(clus.tolist()))

        out.append(sorted(tuple(clusters)))
    return out


def indicatorToCluster(rho):
    out = defaultdict(list)
    for pos, ind in enumerate(rho):
        out[ind].append(pos + 1)

    return json.dumps(sorted(list(out.values())))


def countClusters(rho):
    clusters = indicatorChainToClusters(rho)
    out = defaultdict(int)
    for clus in clusters:
        out[json.dumps(list(clus))] += 1

    return pd.DataFrame([out])


def countNumClusters(rhoChain):
    out = np.zeros(rhoChain.shape[0])
    for i in range(rhoChain.shape[0]):
        out[i] = len(set(rhoChain[i, ].tolist()))

    return out


def areEqualClusters(c1, c2):
    c1 = sorted(tuple([tuple(x) for x in json.loads(c1)]))
    c2 = sorted(tuple([tuple(x) for x in json.loads(c2)]))
    return c1 == c2
