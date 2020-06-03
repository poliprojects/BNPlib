#!/usr/bin/python3
from bnp_interface import *

# Initialize parameters
mu0 = 5.0
lambda_ = 0.1
alpha0 = 2.0
beta0 = 2.0
totalmass = 1.0
algo = "neal2"
init = 3
rng = 20200229
maxit = 5000
burn = 500
n_aux = 3
grid = np.arange(0, 10, 0.1)
only = "all"

# Write file names
datafile = "csv/data_uni.csv"
collfile = "collector.recordio"
densfile = "src/python/dens.csv"
clustfile = "src/python/clust.csv"
imgfileclust = "src/python/clust.pdf"
imgfilechain = "src/python/chain.pdf"
imgfiledens = "src/python/dens.pdf"

# Run algorithms, estimates, and plots
bnplibpy.run_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, datafile, algo,
    collfile, init, rng, maxit, burn, n_aux)

chain_barplot(collfile, imgfilechain)

bnplibpy.estimates_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, grid, algo,
    collfile, densfile, clustfile, only)

plot_clust_cards(clustfile, imgfileclust)
plot_density_points(densfile, imgfiledens)

print("The end")
