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

g = [0.5*_ for _ in range(20)]
grid = np.array(g)

colltype = "file"
collfile = "collector.recordio"
gridfile = "csv/grid_uni.csv"
densfile = "src/python/density.csv"
clustfile = "src/python/best_clust.csv"
imgfilebest = "src/python/best_clust.pdf"
imgfilechain = "src/python/chain.pdf"
imgfiledens = "src/python/density.pdf"

bnplib.run_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, datafile, algo,
	colltype, collfile, rng, maxit, burn)

chain_histogram(collfile, imgfilechain)

algo = "neal8_dataless"
bnplib.estimates_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, gridfile,
	algo, collfile, densfile, clustfile)
# TODO implement estimates cpp s.t. you don't have to write "_dataless"

plot_best_clust_cards(clustfile, imgfilebest)
plot_density(densfile, imgfiledens)
