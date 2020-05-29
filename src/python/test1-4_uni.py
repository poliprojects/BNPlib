#!/usr/bin/python3
from bnp_interface import *

mu0 = 0.0
lambda_ = 0.1
alpha0 = 2.0
beta0 = 2.0
totalmass = 1
algo = "neal2"
algoDL = ''.join((algo, '_dataless'))
rng = 20200229
maxit = 500
burn = 100
colltype = "file"
only = "all"

grid = np.arange(-5, +5.1, 0.5)
#grid2 = np.arange(-7, +7, 0.5)
#grid3 = np.arange(-4, +4, 0.5)
#grid4 = "TODO"
# or linspace

c = [2]
for i in c:
    print("Starting test", i)

    # Write file names
    datafile  = ''.join(("csv/test/data", str(i), ".csv"))
    collfile  = ''.join(("collector", str(i), ".recordio"))
    densfile  = ''.join(("src/python/test_res/dens",  str(i), ".csv"))
    clustfile = ''.join(("src/python/test_res/clust", str(i), ".csv"))
    imgfileclust = ''.join(("src/python/test_res/clust", str(i), ".pdf"))
    imgfilechain = ''.join(("src/python/test_res/chain", str(i), ".pdf"))
    imgfiledens  = ''.join(("src/python/test_res/dens",  str(i), ".pdf"))
    bnplibpy.run_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, datafile,
        algo, colltype, collfile, rng, maxit, burn)

    chain_histogram(collfile, imgfilechain)
    
    bnplibpy.estimates_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, grid,
        algoDL, collfile, densfile, clustfile, only)

    plot_clust_cards(clustfile, imgfileclust)
    plot_density(densfile, imgfiledens)
    trueclustfile = ''.join(("csv/test/true_clust", str(i), ".csv"))
    print_clust_rand_indx(clustfile, trueclustfile)

print("The end")
