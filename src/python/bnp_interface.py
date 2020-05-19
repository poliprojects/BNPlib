#!/usr/bin/python
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import csv
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


def estimates_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, gridfile, algo,
	filecoll_name = "collector.recordio", densfile = "src/python/density.csv"):
	"""TODO docstring

	TODO docstring but longer"""
	bnplib.estimates_NNIG_Dir(mu0, lambda_, alpha0, beta0, totalmass, gridfile,
		algo, filecoll_name, densfile)


def plot_density(densfile = "src/python/density.csv"):
	"""TODO docstring

	TODO docstring but longer"""
	with open(densfile, 'r') as file:
		reader = csv.reader(file, delimiter=',')
		dim = len(next(reader))
		if dim not in (2,3):
			print("Error: density file needs to have either 2 or 3 columns")
			return
		file.seek(0)

		# Get values
		vals = []
		for line in reader:
			vals.append(list(line))
		

		# Plot density
		figure = plt.figure()
		if dim == 2:
			ax = figure.add_subplot(111)
			plt.plot(vals[0], vals[1])
		else: # if dim == 3
			ax = figure.add_subplot(111, projection='3d')
			ax.scatter(vals[0], vals[1], vals[2])
		plt.show()


def histogram():
	pass # plt.hist()
