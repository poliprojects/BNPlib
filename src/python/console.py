#!/usr/bin/python

from bnp_interface import *

mu0 = 5.0
lambda_ = 0.1
alpha0 = 2.0
beta0 = 2.0
totalmass = 1.0
datafile = "csv/data_uni.csv"
algo = "neal2"
coll_type = "memory"
filecoll_name = ""
rng = 20200229
maxit = 1000
burn = 100

bnplib.run_NNIG(mu0, lambda_, alpha0, beta0, totalmass, datafile, algo,
	coll_type, filecoll_name, rng, maxit, burn)
