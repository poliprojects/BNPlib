#!/usr/bin/python
import os, sys
LIBPATH = ''.join((os.path.dirname(os.path.realpath(__file__)), "/../.."))
sys.path.insert(0, LIBPATH)
import bnplib


def run_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, datafile, algo,
	coll_type, filecoll_name = "collector.recordio", rng = 0, maxit = 0,
	burn = 0):
	"""TODO docstring

	TODO docstring but longer"""
	bnplib.run_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, datafile, algo,
		coll_type, filecoll_name, rng, maxit, burn)


def run_NNW_Dir(mu0, lambda_, tau0, nu, totalmass, datafile, algo, coll_type,
	filecoll_name = "collector.recordio", rng = 0, maxit = 0, burn = 0):
	"""TODO docstring

	TODO docstring but longer"""
	bnplib.run_NNW_Dir(mu0, lambda_, tau0, nu, totalmass, datafile, algo,
		coll_type, filecoll_name, rng, maxit, burn)
