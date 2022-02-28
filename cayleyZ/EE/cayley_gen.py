import numpy as np
import scipy.sparse.linalg
import random
import time
import matplotlib.pyplot as plt



def nodes_till_generation(t,z) :
    """
    :param t: generation t
    :param z: coordination number of cayley tree
    :return: nodes till generation t
    """
    if t == 0 :
        res = 0
    else :
        res = int(1 + z*( (z-1)**(t-1) - 1)/(z-2))
    return res


def nodes_in_generation(t,z):
    """
    :param t: genreration t
    :param z: coordination number of cayley tree
    :return:  notes in generation t
    """
    return nodes_till_generation(t,z) - nodes_till_generation(t-1,z)

def find_generation(m,z,total_gen):
    """
    :param m: node of which we need to find generation of
    :param z: coordination number of Cayley tree
    :param total_gen: total generations of the cayley tree
    :return: the generation of the provided node
    """
    for i in range(total_gen) :
        if m <= nodes_till_generation(i+1,z) :
            res = i +1
            break
    return res


def cayley_gen(z,total_gen) :
    """
    :param z: Coordination number of the Cayley tree
    :param total_gen: Generations required of the Cayley tree
    :return: Adjacency matrix of the Cayley tree of the given parameters
    """
    total_nodes = int(1 + z*( (z-1)**(total_gen-1) - 1)/(z-2))
    matrix = np.zeros((total_nodes,total_nodes))

    for m in range(nodes_till_generation(total_gen-1,z)) :
        if m == 0:
            i = 0
            for j in range(1,z+1):
                matrix[i,j] = matrix[i,j] + 1
                matrix[j,i] = matrix[j,i] + 1
        else:
            v = (m+1) - nodes_till_generation(find_generation(m+1,z,total_gen) - 1, z) - 1
            for k in range(1,z):
                i = m
                j = (nodes_till_generation(find_generation(m+1,z,total_gen), z) + (z-1)*v) + k - 1
                matrix[i,j] = matrix[i,j] + 1
                matrix[j,i] = matrix[j,i] + 1
    return matrix

def rvW(W_min,W_max,delta,z,total_gen,cayley_test,loop):
    """
    :param W_min: Minimum strength of disorder
    :param W_max: Maximum strength of disorder
    :param delta: Step size of increasing disorder
    :param z: Coordination number of the Cayley tree
    :param loop: Number of time r is averaged over
    :param total_gen: Number of generations in the Cayley tree, central vertex is taken as generation 1.
    :param cayley_test: Adjacency matrix of the Cayley tree
    :return: Array of W, r (Ensemble averaged ratio) and error (standard deviation scaled as 1/10 of original)
    """
    total_vert = nodes_till_generation(total_gen,z)     #Total vertices in the Cayley tree
    avgrrange = np.array([])
    avgavgrrange = np.array([])
    err = np.array([])
    for k in range(loop):
        #print("ok - " + str(k))
        #print("ok1 - " + str(k))
        W = W_min
        # print(np.linalg.eig(H))
        for j in range(int(W_max/delta)):
            W = W + delta
            for i in range(total_vert) :
                w = W*(random.random() - 0.5)
                cayley_test[i,i] = w


            #eigval = scipy.sparse.linalg.eigs(H, k = int(n/2), sigma=0, return_eigenvectors=False)
            eigval,eigvec = scipy.linalg.eigh(cayley_test,overwrite_a=False,check_finite=True,turbo=True)
            sorted_eigval = np.sort(eigval)
            level_spacing = np.array([])
            for i in range(total_vert-1):
                level_spacing = np.append(level_spacing,sorted_eigval[i+1]-sorted_eigval[i])
            #print(level_spacing)

            rn = np.array([])

            for i in range(total_vert-2) :
                rn = np.append(rn,min(level_spacing[i],level_spacing[i+1])/max(level_spacing[i],level_spacing[i+1]))

            avgr = sum(rn)/len(rn)
            #print(avgr)
            avgrrange = np.append(avgrrange,avgr)

        if k == 0 :
            avgavgrrange = avgrrange.copy()
        else :
            avgavgrrange = np.vstack([avgavgrrange,avgrrange])
        avgrrange = np.array([])

    for i in range(len(avgavgrrange[0,:])):
        err = np.append(err,np.std(avgavgrrange[:,i])/10)
    finalr = np.mean(avgavgrrange,axis=0)
    W = np.linspace(W_min,W_max,int((W_max-W_min)/delta))
    return W,finalr,err

def wavefunc_t(cayley,total_vert,initial,t):
    eigval,eigvec = scipy.linalg.eigh(cayley,overwrite_a=True,check_finite=True,turbo=True)
    #eigval,eigvec = np.linalg.eig(H)
    #eigval,eigvec = scipy.sparse.linalg.eigsh(H, k = int(n/2), sigma=0, return_eigenvectors=True)
    psi = np.zeros(total_vert)
    for i in range(total_vert):
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

def ipr(psi):
    """
    :param psi: wavefunction psi
    :return: ipr
    """
    #final = wavefunc_t(cayley_cpy,total_vert,initial,time_val)
    #final = abs(final)
    psi = abs(psi)
    for i in range(len(psi)):
        psi[i] = psi[i]**4
    return sum(psi)
def neighbours_in_nearest_higher_generation_till_end(a,cayley):
    """
    :param a: initial vertex as array (center vertex is one)
    :return: neighbours_in_nearest_higher_generation_till_end - a
    """
    a = a - 1
    b = np.array([])
    condition = False
    c = np.array([3])
    while condition == False:
        while True:
            for i in a:
                b = np.append(b,np.array(np.nonzero(cayley[int(i),int(i+1):])) + int(i+1))
            b = np.unique(b)
            if len(c) == len(b):
                break
            c = b
            for item in b:
                if item not in a:
                    a = np.append(a,item)
        condition = True
    return a + 1

def nodes_in_gen_bond_center(g):
    """
    :param g: generation in coordination number 3
    :return:
    """
    return int(2**g)

def nodes_till_gen_bond_center(g):
    """
    :param g: generation (coordination number is 3)
    :return: nodes till generation g including g
    """
    count = 0
    for i in range(g):
        count = 2**(i+1) + count
    return int(count)

def cayley_gen_bond_center(G):
    t =nodes_till_gen_bond_center(G)
    matrix = np.zeros((t,t))
    t = nodes_till_gen_bond_center(G-1)
    matrix[0,2**G - 1] = 1
    matrix[2 ** G - 1,0] = 1
    for i in range(int(t/2)):
        matrix[i , 2 * ( i + 1 ) -1] = 1
        matrix[i, 2 * (i + 1)] = 1
        matrix[2 * (i + 1) - 1, i] = 1
        matrix[2 * (i + 1), i] = 1
        matrix[2**G - 1 + i , 2**G - 1 + 2 * ( i + 1 ) -1] = 1
        matrix[2**G - 1 + i, 2**G - 1 + 2 * (i + 1)] = 1
        matrix[2**G - 1 + 2 * (i + 1) - 1, 2**G - 1 + i] = 1
        matrix[2**G - 1 + 2 * (i + 1), 2**G - 1 + i] = 1
    return matrix