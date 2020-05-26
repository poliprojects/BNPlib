#!/usr/bin/python
from bnp_interface import *
from multi_utils import *

lambda_ = 0.2;
totalmass = 1
algo = "neal2"
algoDL = ''.join((algo, '_dataless'))
rng = 20200229
maxit = 1000
burn = 100
colltype = "file"
only = "all"


d=[2,5,10,20]
c=[5] # TODO per plot
for i in c:
	print("Starting test", i)
	datafile  = ''.join(("csv/test/data", str(i), ".csv"))
	collfile  = ''.join(("collector", str(i), ".recordio"))
	densfile  = ''.join(("src/python/test_res/dens",  str(i), ".csv"))
	clustfile = ''.join(("src/python/test_res/clust", str(i), ".csv"))
	imgfileclust = ''.join(("src/python/test_res/clust", str(i), ".pdf"))
	imgfilechain = ''.join(("src/python/test_res/chain", str(i), ".pdf"))
	imgfiledens  = ''.join(("src/python/test_res/dens",  str(i), ".pdf"))
	mu0 = empirical_mean(datafile)
	nu=d[i-5]+3
	tau0 = (1/nu) * np.identity(d[i-5])
	bnplibpy.run_NNW_Dir(mu0, lambda_,tau0, nu, totalmass, datafile, algo,colltype, collfile, rng, maxit, burn)
	chain_histogram(collfile, imgfilechain)
	grid=get_grid(d[i-5])
	bnplibpy.estimates_NNW_Dir(mu0, lambda_, tau0, nu, totalmass, grid[0], algoDL,
	collfile, densfile, clustfile, only)
	plot_clust_cards(clustfile, imgfileclust)
	plot_density(densfile, imgfiledens)




