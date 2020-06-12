##
# @file
# Tools for interfacing the C++ BNPlib library and Python
#

from google.protobuf.internal.encoder import _VarintBytes
from google.protobuf.internal.decoder import _DecodeVarint32
import chain_state_pb2

import csv, os, sys
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from sklearn.metrics.cluster import adjusted_rand_score

LIBPATH = ''.join((os.path.dirname(os.path.realpath(__file__)), "/../.."))
sys.path.insert(0, LIBPATH)
import bnplibpy


def get_multidim_grid(a, b, d, n=10):
    """! Builds a d-dimensional grid from intervals [a,b], each divided into n.

    Given the extrema a and b, this function creates a one-dimensional array
    with n sample points, then builds coordinate matrices from the single d
    coordinate vectors, and returns a hypercubic grid where every matrix is
    collapsed into one dimension."""
    uni_g = np.linspace(a, b, n)
    arr = [uni_g for y in range(d)]
    mesh = np.meshgrid(*arr)
    nrows = len(mesh[0].flat)
    grid = np.zeros((nrows, d))
    for i in range(0,d):
        grid[:,i] = mesh[i].flat
    return grid


def deserialize(collfile):
    """! Reads collector from file and returns it as a list of Protobuf objects

    collfile is the name of the saved file collector. All Protobuf messages in
    the returned list are instances of the custom State class, that is, they
    correspond to the whole state of a single iteration of an algorithm of this
    library."""
    with open(collfile, 'rb') as f:
        buf = f.read()
        n = 0
        d = []
        while n < len(buf):
            msg_len, new_pos = _DecodeVarint32(buf, n)
            n = new_pos
            msg_buf = buf[n:n+msg_len]
            n += msg_len
            read_metric = chain_state_pb2.State()
            read_metric.ParseFromString(msg_buf)
            d.append(read_metric)
    return d


def chain_barplot(collfile, imgfile = "src/python/chain.pdf"):
    """! Prints a barplot of the number of clusters for different iterations.

    collfile is the name of the saved file collector. After deserialization, for
    each State object (i.e. for every iteration of the algorithm) the number of
    clusters at that time is counted and plotted in a barplot, which is then
    saved to the imgfile file."""
    d = deserialize(collfile)
    figure = plt.figure()
    num_clusters = []
    for i in d:
        num_clusters.append(len(i.uniquevalues))

    hist, bin_edges = np.histogram( num_clusters,
        bins=range(min(num_clusters), max(num_clusters)+2) )
    plt.bar(bin_edges[:-1], hist)
    plt.title("Frequency of number of clusters over iterations")
    plt.xlabel("Number of clusters")
    plt.ylabel("Absolute frequency")
    plt.xticks(bin_edges[:-1])
    plt.savefig(imgfile)
    print("Successfully saved plot to", imgfile)


def plot_density_points(densfile, imgfile = "src/python/dens_points.pdf",
    title = "Posterior estimate graph"):
    """! Reads a 1D or 2D density from a csv file and plots it point-by-point.

    densfile is the file from which the density will be read. Such file must
    consist of 2 to 3 columns, the last of which is the density value and the
    remaining ones are the coordinates of the 1D or 2D points of evaluation of
    that value. Then, a 2D or 3D plot with the given title is produced and saved
    to the imgfile file. If instead the columns are neither 2 nor 3, a plot is
    not produced."""
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
    plt.title(title)
    plt.savefig(imgfile)
    print("Successfully saved plot to", imgfile)


def plot_density_contour(densfile, imgfile = "src/python/dens_cont.pdf",
    title = "Posterior estimate contour"):
    """! Reads a 2D density from a csv file and plots countour lines for it.

    densfile is the file from which the density will be read. Such file must
    consist of 3 columns, the last of which is the density value while the first
    2 are the coordinates of the 2D points of evaluation of that value. Then, a
    2D plot with the given title is produced and saved to the imgfile file. If
    instead the columns are not 3, a plot is not produced."""
    mat = np.loadtxt(open(densfile, 'rb'), delimiter=',')
    cols = mat.shape[1]
    if cols != 3:
        print("Error: density file must have 3 columns to be plotted")
        return
    figure = plt.figure()
    n = np.floor(np.sqrt(mat[:,0].shape[0]))
    x = np.linspace(min(mat[:,0]), max(mat[:,0]), n.astype(int))
    xx, yy = np.meshgrid(x, x)
    z = mat[:,2].reshape(xx.shape)
    plt.contourf(xx, yy, z)
    plt.title(title)
    plt.savefig(imgfile)
    print("Successfully saved plot to", imgfile)


def plot_clust_cards(clustfile, imgfile = "src/python/clust.pdf"):
    """! Reads a data clustering from a csv file and plots the cardinalities.

    clustfile is the file from which the clustering will be read. This function
    reads the first column of the file and interprets it as the numeric labels
    of the clusters the data points belong to, with each for being a different
    datum. A barplot of the cardinalities of these clusters is then produced and
    saved to the imgfile file."""
    mat = np.loadtxt(open(clustfile, 'rb'), delimiter=',')
    figure = plt.figure()
    clusters = mat[:,0].astype(int)
    hist, bin_edges = np.histogram( clusters,
        bins=range(min(clusters), max(clusters)+2) )
    plt.bar(bin_edges[:-1], hist)
    plt.xticks(bin_edges[:-1])
    plt.title("Cardinalities of clustering")
    plt.xlabel("Cluster label")
    plt.ylabel("Cardinality")
    plt.savefig(imgfile)
    print("Successfully saved plot to", imgfile)


def clust_rand_score(clustfile, trueclustfile):
    """! Returns the adjusted Rand score between predicted and true clustering.

    clustfile and trueclustfile are the file names from which the predicted and
    the true clusterings, respectively, will be read. The cluster labels must be
    in the first column of the csv files. Then, this function computes and
    returns the adjusted Rand score of the predicted clustering, which is a
    measure of the proportion of correct decisions made by the clustering
    algorithm."""
    mat_pred = np.loadtxt(open(clustfile, 'rb'), delimiter=',')
    mat_true = np.loadtxt(open(trueclustfile, 'rb'), delimiter=',')
    score = adjusted_rand_score(mat_pred[:,0].astype(int), mat_true.astype(int))
    return score
