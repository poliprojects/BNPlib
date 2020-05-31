#!/usr/bin/python3
from bnp_interface import *

lambda_ = 0.2
totalmass = 1
algo = "neal2"
rng = 20200229
maxit = 200
burn = 10
n_aux = 3
only = "all"

d = [2,5,10,20]
c = [5] # TODO per plot
for i in c:
    print("Starting test", i)
    datafile  = ''.join(("csv/test/data", str(i), ".csv"))
    collfile  = ''.join(("collector", str(i), ".recordio"))
    densfile  = ''.join(("src/python/test_res/dens",  str(i), ".csv"))
    clustfile = ''.join(("src/python/test_res/clust", str(i), ".csv"))
    imgfileclust = ''.join(("src/python/test_res/clust", str(i), ".pdf"))
    imgfilechain = ''.join(("src/python/test_res/chain", str(i), ".pdf"))
    imgfiledens  = ''.join(("src/python/test_res/dens",  str(i), ".pdf"))
    imgfilecontour  = ''.join(("src/python/test_res/dens_contour",  str(i),
        ".pdf"))
    mat = np.loadtxt(open(datafile, 'rb'), delimiter=' ')
    mu0 = np.mean(mat, axis=0)
    nu = d[i-5]+3
    tau0 = (1/nu) * np.identity(d[i-5])
    bnplibpy.run_NNW_Dir(mu0, lambda_,tau0, nu, totalmass, datafile, algo,
        collfile, rng, maxit, burn, n_aux)
    chain_histogram(collfile, imgfilechain)
    grid = get_grid(d[i-5])
    bnplibpy.estimates_NNW_Dir(mu0, lambda_, tau0, nu, totalmass, grid, algo,
    	collfile, densfile, clustfile, only)
    plot_clust_cards(clustfile, imgfileclust)
    plot_density_points(densfile, imgfiledens)
    plot_density_contour(densfile, imgfilecontour)
