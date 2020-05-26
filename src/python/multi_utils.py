#!/usr/bin/python
import numpy as np

def get_grid(d):
	uni_g = np.arange(-5, +5.1, 0.5)
	arr = [uni_g for y in range(d)]
	mesh= np.meshgrid(*arr)
	nrows=len(mesh[0].flat)
	grid = np.zeros((nrows,d))


	for i in range(0,d):
		grid[:,i] = mesh[i].flat
		
	return[grid]


def empirical_mean(datafile):
	mat = np.loadtxt(open(datafile, 'rb'), delimiter=' ')
	return[np.mean(mat, axis=0)]

