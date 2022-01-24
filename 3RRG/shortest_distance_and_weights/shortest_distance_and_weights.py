import numpy as np
from Making_graph_matrix import *
import math
import scipy
import scipy.linalg
import matplotlib.pyplot as plt


def nearest_neighbour(initial,rrg):
    """
    :param initial: set of initial vertices of whose the neighbours needs to be found. Indexing starts from 0
    :param rrg: adjacency matrix of the random regular graph.
    :return: the neighbouring vertices of the vertices provided in arg.

    """
    neighbour_set = np.array([])
    for j in range(len(initial)):
        for i in range(n):
            if rrg[(int(initial[j])-1),i] == 1:
                neighbour_set = np.append(neighbour_set,i+1)

    return neighbour_set

def distance_finder(initial,end,rrg):
    """
    :param initial: the numpy array containing initial vertex(also can be a set of vertices)
    :param end: the vertex index to which the distance from the initial vertex needs to be found
    :param rrg: adjacency matrix of the random regular graph
    :return: the minimum distance/number of edges from the initial vertex to the end vertex
    """
    #print("distance_finder has started")
    d = 0
    condition = False
    if (end in initial) == True:
        condition = True
    while condition == False:
        if d == 0:
            neighbours = nearest_neighbour(initial,rrg)
            d += 1
        elif (end in neighbours) == True:
            condition = True
        else:
            neighbours = nearest_neighbour(neighbours,rrg)
            d += 1
    #print("distance_finder has stopped")
    return d



def wavefunc_t(rrg,n,initial,t):
    eigval,eigvec = scipy.linalg.eigh(rrg,overwrite_a=True,check_finite=True,turbo=True)
    #eigval,eigvec = np.linalg.eig(H)
    #eigval,eigvec = scipy.sparse.linalg.eigsh(H, k = int(n/2), sigma=0, return_eigenvectors=True)
    psi = np.zeros(n)
    for i in range(n):
        psi = psi + math.e**(-1j*t*eigval[i])*(np.dot(np.transpose(eigvec[:,i]),initial))*eigvec[:,i]
        #psi = psi + sp.exp((-1j*t*eigval[i]))*((sp.transpose(sp.Matrix(eigvec[:,i])).dot(sp.Matrix(initial))))*sp.Matrix(eigvec[:,i])
    return psi

def prob_density(psi):
    #final = wavefunc_t(cayley_cpy,total_vert,initial,time_val)
    #final = abs(final)
    psi = abs(psi)
    for i in range(len(psi)):
        psi[i] = psi[i]**2
    return psi

def hist_weight(psi,n,initial_index,rrg):
        """
        :param psi: wavefuntion amplitudes
        :param n: total number of nodes
        :param initial: initial vertex index of particle in numpy array
        :param rrg: the random regular graph in question
        :return: weight of the graph at given time wavefunction
        """
        #print("weight calculation has started")
        psi = prob_density(psi)
        hist = np.zeros((n,1))

        for i in range(n):
            initial_cpy = initial_index.copy()
            hist[i] = distance_finder(initial_cpy,i+1,rrg)*psi[i]
        hist_weight_sum = sum(hist)
        #print("weight calculation has stopped")

        return hist_weight_sum

def weightFuncOf_t_vs_Time(n,rrg,initial,initial_index,total_time,W):
    # weight is the summation of probability magnitude at the lattice site times the distance of that lattice from the
    # initial vertex. The initial vertex is given zero weight, i.e. more the spread, more the weight.
    for i in range(n) :
        w = W*(np.random.random() - 0.5)
        rrg[i,i] = w
    weight_func_t = np.array([])
    for i in range(total_time):
        weight_func_t = np.append(weight_func_t,hist_weight(wavefunc_t(rrg,n,initial,i/10),n,initial_index,rrg))
    time_points = np.linspace(0,total_time-1,total_time)
    time_points = time_points/10
    return weight_func_t,time_points


############################
# Multiple subplots for given
# disorder over a range
############################
time_val = 101
W = 0
n = 128
initial_index = np.array([32])
initial = np.zeros((n,1))
initial[31] = 1
rrg = three_rrg(n)
plt.rcParams['axes.grid'] = True
for k in range(15):
    W = W + 2
    for j in range(5):
        plt.title(f'Weight plots; W = {W}; N = {n}')
        weight_func_t,time_points = weightFuncOf_t_vs_Time(n,rrg,initial,initial_index,time_val,W)
        plt.plot(time_points,weight_func_t,'r-')
        plt.xlabel("Time")
        plt.ylabel("Weight")
        plt.ylim(0,5)
        plt.tight_layout()
        plt.savefig(f"weight_func_vs_t_W={W}_({j+1})_N={n}.png", format="png", dpi=1280)
        #plt.show()
        plt.clf()
