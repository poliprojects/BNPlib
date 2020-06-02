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
from scipy.special import comb

LIBPATH = ''.join((os.path.dirname(os.path.realpath(__file__)), "/../.."))
sys.path.insert(0, LIBPATH)
import bnplibpy


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


def chain_histogram(collfile, imgfile = "src/python/chain.pdf"):
    """! Prints an histogram of the number of clusters for different iterations.

    collfile is the name of the saved file collector. After deserialization, for
    each State object (i.e. for every iteration of the algorithm) the number of
    clusters at that time is counted and plotted in an histogram, which is then
    saved to the imgfile file."""
    d = deserialize(collfile)
    figure = plt.figure()
    num_clusters = []
    for i in d:
        num_clusters.append(len(i.uniquevalues))
    plt.hist(num_clusters)
    plt.savefig(imgfile)
    print("Successfully saved plot to", imgfile)


def plot_density_points(densfile, imgfile = "src/python/dens_points.pdf"):
    """! Reads a 1D or 2D density from a csv file and plots it point-by-point.

    densfile is the file from which the density will be read. Such file must
    consist of 2 to 3 columns, the last of which is the density value and the
    remaining ones are the coordinates of the 1D or 2D points of evaluation of
    that value. Then, a 2D or 3D plot is produced and saved to the imgfile file.
    If instead the columns are neither 2 nor 3, a plot is not produced."""
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


def plot_density_contour(densfile, imgfile = "src/python/dens_cont.pdf"):
    """! Reads a 2D density from a csv file and plots countour lines for it.

    densfile is the file from which the density will be read. Such file must
    consist of 3 columns, the last of which is the density value while the first
    2 are the coordinates of the 2D points of evaluation of that value. Then, a
    2D plot is produced and saved to the imgfile file. If instead the columns
    are not 3, a plot is not produced."""
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
    plt.savefig(imgfile)
    print("Successfully saved plot to", imgfile)


def plot_clust_cards(clustfile, imgfile = "src/python/clust.pdf"):
    """! Reads a data clustering from a csv file and plots the cardinalities.

    clustfile is the file from which the clustering will be read. This function
    reads the first column of the file and interprets it as the numeric labels
    of the clusters the data points belong to, with each for being a different
    datum. An histogram of the cardinalities of these clusters is then produced
    and saved to the imgfile file."""
    mat = np.loadtxt(open(clustfile, 'rb'), delimiter=',')
    figure = plt.figure()
    plt.hist(mat[:,0])
    plt.savefig(imgfile)
    print("Successfully saved plot to", imgfile)


def rand_index_score(clusters, classes):
    """! Computes the Rand index score between predicted and true clustering.

    The inputs are 2 lists of labels; clusters, the predicted clustering, and
    classes, the true clustering. After counting all true/false positives/
    negatives, the Rand index is returned, which is a measure of the percentage
    of correct decisions made by the algorithm, computed as the empirical
    proportion of the correctly classified points."""
    tp_plus_fp = comb(np.bincount(clusters), 2).sum()
    tp_plus_fn = comb(np.bincount(classes), 2).sum()
    A = np.c_[(clusters, classes)]
    tp = sum( comb(np.bincount(A[A[:, 0] == i, 1]), 2).sum()
        for i in set(clusters) )
    fp = tp_plus_fp - tp
    fn = tp_plus_fn - tp
    tn = comb(len(A), 2) - tp - fp - fn
    return (tp + tn) / (tp + fp + fn + tn)


def print_clust_rand_index(clustfile, trueclustfile):
    """! Computes the Rand index score between predicted and true clustering.

    clustfile and trueclustfile are the file names from which the predicted and
    the true clusterings, respectively, will be read. The cluster labels must be
    in the first column of the csv files, and clustfile must also contain the
    unique values that characterize the clusters. This function finds the k most
    numerous predicted clusters, where k is the number of true clusters, and
    relabels them. Relabeling is done in ascending order with respect to the
    first unique value (which is the cluster mean in most cases), assuming that
    the true clustering file is sorted in the same way. Then, the Rand index
    score is computed and printed using the function. Please note that this
    index is usually computed when the number of clusters is known a priori and
    fixed, and choosing only the biggest k clusters as significant (while
    misclassifying all other points) is an approximation; a useful approximation
    nonetheless."""
    mat_pred = np.loadtxt(open(clustfile, 'rb'), delimiter=',')
    mat_true = np.loadtxt(open(trueclustfile, 'rb'), delimiter=',')

    clusters, tot_counts = np.unique(mat_pred[:,0].astype(int),
        return_counts=True)
    biggest_counts = sorted(tot_counts, reverse=True)[0:len(np.unique(
        mat_true))]
    biggest_clust = clusters[ np.where(np.in1d(tot_counts, biggest_counts))[0]]

    v1 = np.unique(mat_pred[np.where(np.in1d(mat_pred[:,0],biggest_clust))[0],
        1])
    v2 = np.array([item for item in np.unique(mat_pred[:,1]) if item not in v1])
    values = np.concatenate((v1,v2))

    new_indx = np.arange(len(values))

    for i in range(0, mat_pred.shape[0]):
        mat_pred[i,0] = new_indx[ np.where(values == mat_pred[i,1]) ]

    print("Rand index score:", rand_index_score(mat_pred[:,0].astype(int),
        mat_true.astype(int)))
