#!/usr/bin/python3
from bnp_interface import *

# Initialize parameters
lambda_ = 0.2
totalmass = 1
algo = "neal2"
init = 0
rng = 20200229
maxit = 200
burn = 10
n_aux = 3
only = "all"

dim = [2,5]
#dim = [2,5,10,20]
tests = [5,6]

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
    imgfilecont   = ''.join(("src/python/test_res/cont",  str(t), ".pdf"))

    # Initialize more parameters
    mat = np.loadtxt(open(datafile, 'rb'), delimiter=' ')
    mu0 = np.mean(mat, axis=0)
    nu = dim[t-5] + 3
    tau0 = (1/nu) * np.identity(dim[t-5])

    bnplibpy.run_NNW_Dir(mu0, lambda_,tau0, nu, totalmass, datafile, algo,
        collfile, init, rng, maxit, burn, n_aux)

    chain_histogram(collfile, imgfilechain)

    grid = get_multidim_grid(-7, 7.1, dim[t-5], 10)

    bnplibpy.estimates_NNW_Dir(mu0, lambda_, tau0, nu, totalmass, grid, algo,
    	collfile, densfile, clustfile, only)

    plot_clust_cards(clustfile, imgfileclust)
    plot_density_points(densfile, imgfiledens)
    plot_density_contour(densfile, imgfilecont)

    print()
