#!/usr/bin/python
print("Importing...")
from bnp_interface import *

print("Done")
mu0 = 5.0
lambda_ = 0.1
alpha0 = 2.0
beta0 = 2.0
totalmass = 1.0
datafile = "csv/data_uni.csv"
#grid = [0.5*_ for _ in range(20)]
gridfile = "csv/grid_uni.csv"



algo = "neal8"
coll_type = "file"
filecoll_name = "collector.recordio"
densfile = "src/python/density.csv"
rng = 20200229
maxit = 1000
burn = 100

bnplib.run_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, datafile, algo,
	coll_type, filecoll_name, rng, maxit, burn)

algo = "neal8_dataless"
bnplib.estimates_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, gridfile,
	algo, filecoll_name, densfile)

plot_density(densfile)

#plot_density("csv/dens_multi.csv")
