#!/usr/bin/python3
from bnp_interface import *

# Initialize parameters
mu0 = 0
lambda_ = 0.1
alpha0 = 2.0
beta0 = 2.0
totalmass = 1
algo = "neal2"
init = 0
rng = 20200229
maxit = 500
burn = 100
n_aux = 3
only = "all"

# Build grid for density evaluation
grids = []
grids.append( np.arange( -7,  +7.1, 0.1) ) # for test 1
grids.append( np.arange(-10, +10.1, 0.1) ) # for test 2
grids.append( np.arange( -6,  +6.1, 0.1) ) # for test 3
grids.append( np.arange(-10, +10.1, 0.1) ) # for test 4

tests = [1,2,3,4]

for t in tests:
    print("Starting test", t)

    # Write file names
    datafile  = ''.join(("csv/test/data", str(t), ".csv"))
    collfile  = ''.join(("collector", str(t), ".recordio"))
    densfile      = ''.join(("src/python/test_res/dens",  str(t), ".csv"))
    clustfile     = ''.join(("src/python/test_res/clust", str(t), ".csv"))
    imgfileclust  = ''.join(("src/python/test_res/clust", str(t), ".pdf"))
    imgfilechain  = ''.join(("src/python/test_res/chain", str(t), ".pdf"))
    imgfiledens   = ''.join(("src/python/test_res/dens",  str(t), ".pdf"))
    trueclustfile = ''.join(("csv/test/true_clust", str(t), ".csv"))

    # Run algorithms, estimates, and plots
    bnplibpy.run_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, datafile,
        algo, collfile, init, rng, maxit, burn, n_aux)

    chain_barplot(collfile, imgfilechain)

    bnplibpy.estimates_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass,
        grids[t-1], algo, collfile, densfile, clustfile, only)

    plot_clust_cards(clustfile, imgfileclust)
    plot_density_points(densfile, imgfiledens)
    print("Adjusted Rand score:", clust_rand_score(clustfile,trueclustfile))

    print()
