from google.protobuf.internal.encoder import _VarintBytes
from google.protobuf.internal.decoder import _DecodeVarint32
import output_pb2

import os, sys
import matplotlib.pyplot as plt
import numpy as np
import csv
from mpl_toolkits.mplot3d import Axes3D
from sklearn.metrics.cluster import adjusted_rand_score

LIBPATH = ''.join((os.path.dirname(os.path.realpath(__file__)), "/../.."))
sys.path.insert(0, LIBPATH)
import bnplibpy


def deserialize(collfile = "collector.recordio"):
	with open(collfile, 'rb') as f:
		buf = f.read()
		n = 0
		d = []
		while n < len(buf):
			msg_len, new_pos = _DecodeVarint32(buf, n)
			n = new_pos
			msg_buf = buf[n:n+msg_len]
			n += msg_len
			read_metric = output_pb2.State()
			read_metric.ParseFromString(msg_buf)
			d.append(read_metric)
	return[d]


def chain_histogram(collfile = "collector.recordio",
	imgfile = "src/python/chain.pdf"):
	d = deserialize(collfile)
	figure = plt.figure()
	num_clusters = []
	
	for i in d[0]:  
		num_clusters.append(len(i.uniquevalues))

	plt.hist(num_clusters)
	plt.savefig(imgfile)
	print("Successfully saved plot to", imgfile)


def plot_density(densfile = "src/python/density.csv",
	imgfile = "src/python/density.pdf"):
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
	print("Successfully saved plot to", imgfile)


def plot_clust_cards(clustfile = "src/python/clust.csv",
	imgfile = "src/python/clust.pdf"):
	mat = np.loadtxt(open(clustfile, 'rb'), delimiter=',')
	figure = plt.figure()
	plt.hist(mat[:,0])
	plt.savefig(imgfile)
	print("Successfully saved plot to", imgfile)
	
	
def print_clust_rand_indx(clustfile = "src/python/clust.csv",
	trueclustfile = "src/csv/test/true_clust.csv"):
	mat = np.loadtxt(open(clustfile, 'rb'), delimiter=',')
	mat_true = np.loadtxt(open(trueclustfile, 'rb'), delimiter=',')
	print(adjusted_rand_score(mat[:,0], mat_true))
