#!/usr/bin/python
from bnp_interface import *

mu0 = 5.0
lambda_ = 0.1
alpha0 = 2.0
beta0 = 2.0
totalmass = 1.0
datafile = "csv/data_uni.csv"
algo = "neal8"
rng = 20200229
maxit = 5000
burn = 500

coll_type = "file"
filecoll_name = "collector.recordio"

g= [0.5*_ for _ in range(20)]
grid = np.array(g)

gridfile = "csv/grid_uni.csv"
densfile = "src/python/density.csv"
imgfile  = "src/python/plot.pdf"

bnplib.run_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, datafile, algo,
	coll_type, filecoll_name, rng, maxit, burn)

algo = "neal8_dataless"
bnplib.estimates_NNIG_Dir_grid(mu0, lambda_, alpha0, beta0, totalmass, grid,
	algo, filecoll_name, densfile)
# TODO implement estimates cpp s.t. you don't have to write "_dataless"

plot_density(densfile, imgfile)
