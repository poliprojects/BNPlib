#!/usr/bin/python3
from bnp_interface import *

mu0 = 5.0
lambda_ = 0.1
alpha0 = 2.0
beta0 = 2.0
totalmass = 1.0
datafile = "csv/data_uni.csv"
algo = "neal2"
rng = 20200229
maxit = 5000
burn = 500

g = [0.5*_ for _ in range(20)]
grid = np.array(g)

colltype = "file"
collfile = "collector.recordio"
densfile = "src/python/dens.csv"
clustfile = "src/python/clust.csv"
imgfileclust = "src/python/clust.pdf"
imgfilechain = "src/python/chain.pdf"
imgfiledens = "src/python/dens.pdf"
only = "all"

bnplibpy.run_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, datafile, algo,
    colltype, collfile, rng, maxit, burn)

chain_histogram(collfile, imgfilechain)

algo = ''.join((algo, '_dataless'))
bnplibpy.estimates_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, grid, algo,
    collfile, densfile, clustfile, only)
# TODO implement estimates cpp s.t. you don't have to write "_dataless"
# TODO: is it really worth keeping the dataless constructors in the factory?

plot_clust_cards(clustfile, imgfileclust)
plot_density(densfile, imgfiledens)
g = np.arange(-5, +5, 0.5)
grid = np.repeat(g,2).T
matr = np.linspace(np.repeat(0, 2),np.repeat(10, 2),10)

z = np.identity(10)
print(z)

print("The end")
