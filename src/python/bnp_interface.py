"""@package docstring
BNP Python interface.
 
Documentation for BNP Python interface, but longer.
"""
from deserialize import deserialize
import output_pb2

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


def estimates_NNIG_Dir_grid(mu0, lambda_, alpha0, beta0, totalmass, grid, algo,
	filecoll_name = "collector.recordio", densfile = "src/python/density.csv"):
	"""TODO docstring

	TODO docstring but longer"""
	bnplib.estimates_NNIG_Dir_grid(mu0, lambda_, alpha0, beta0, totalmass, grid,
		algo, filecoll_name, densfile)


def plot_density(densfile = "src/python/density.csv",
	imgfile = "src/python/plot.pdf"):
	"""TODO docstring.

	TODO docstring but longer."""
	mat = np.loadtxt(open(densfile, 'rb'), delimiter=',')


	cols = mat.shape[1]
	if cols not in (2,3):
		print("Error: density file must have 2-3 columns to be plotted")
		return
	figure = plt.figure()
	if cols == 2:
		ax = figure.add_subplot(111)
		plt.plot(mat[:,0], mat[:,1])
	else: # if cols == 3
		ax = figure.add_subplot(111, projection='3d')
		ax.scatter(mat[:,0], mat[:,1], mat[:,2])
	plt.savefig(imgfile)


def histogram(imgfile = "src/python/hist.pdf"):
	
	d=deserialize()
	figure = plt.figure()
	num_clusters=[]
	
	for i in d[0]:  
  		num_clusters.append(len(i.uniquevalues))

	plt.hist(num_clusters)
	plt.savefig(imgfile)	
	
	
	
def plotcardinalities(clusterfile= "src/python/best_clustering.csv", imgfile = "src/python/best_clustering.pdf"):
	mat = np.loadtxt(open(clusterfile, 'rb'), delimiter=',')
	figure = plt.figure()
	plt.hist(mat[:,0])
	plt.savefig(imgfile)	
	
	
	
	
