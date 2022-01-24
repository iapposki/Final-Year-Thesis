from cayley_gen import *
import numpy as np
import math
import scipy.sparse.linalg
import random
import time
import matplotlib.pyplot as plt
import sympy as sp


begin = time.time()

total_gen = 4
z = 3
W = 20
time_val = sp.symbols('time_val')

total_vert = nodes_till_generation(total_gen,z)
initial = np.zeros((total_vert,1))
initial[0] = 1
cayley = cayley_gen(z,total_gen)
for i in range(total_vert) :
    w = W*(random.random() - 0.5)
    cayley[i,i] = w
cayley_cpy = cayley.copy()
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

def hist(psi,total_gen,z):
    """
    :param psi: wavefunction amplitudes.
    :param total_gen: total generations.
    :param z: coordination number.
    :return: probability density apmlitudes as function of generations.
    """
    psi = prob_density(psi)
    hist = np.zeros((total_gen,1))

    for i in range(total_gen):
        if i == 0:
            hist[0] = psi[0]
        else:
            hist[i] = sum(psi[nodes_till_generation(i-1,z)+1:nodes_till_generation(i,z)+1])
    return hist



n_pd = np.linspace(1,total_vert,total_vert)
n_hist = np.linspace(1,total_gen,total_gen)


which = input('Which(pd/hist)?')

if which == 'pd':
    for i in range(210):
        plt.plot(n_pd,prob_density(wavefunc_t(cayley_cpy,total_vert,initial,i/10)),'r-')
        plt.xticks(np.arange(min(n_pd), max(n_pd)+1, 1.0))
        plt.title(f'time = {i/10} ')
        plt.ylim(0,1)
        plt.grid()
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        plt.show()
elif which == 'hist':
    for i in range(210):
        plt.plot(n_hist,hist(wavefunc_t(cayley_cpy,total_vert,initial,i/10),total_gen,z),'r-')
        plt.grid()
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        plt.show()

end = time.time()
print("Time taken is " + str(end-begin) + " seconds")










