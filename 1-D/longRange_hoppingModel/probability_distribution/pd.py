import numpy as np
import matplotlib.pyplot as plt
import math
import random
import scipy.linalg



def wavefunc_t(matrix,n,initial,t):
    eigval,eigvec = scipy.linalg.eigh(matrix,overwrite_a=True,check_finite=True,turbo=True)
    #eigval,eigvec = np.linalg.eig(H)
    #eigval,eigvec = scipy.sparse.linalg.eigsh(H, k = int(n/2), sigma=0, return_eigenvectors=True)
    psi = np.zeros(n)
    for i in range(n):
        psi = psi + math.e**(-1j*t*eigval[i])*(np.dot(np.transpose(eigvec[:,i]),initial))*eigvec[:,i]
        #psi = psi + sp.exp((-1j*t*eigval[i]))*((sp.transpose(sp.Matrix(eigvec[:,i])).dot(sp.Matrix(initial))))*sp.Matrix(eigvec[:,i])
    return psi

def prob_density(psi):
    #final = wavefunc_t(matrix_cpy,n,initial,time_val)
    #final = abs(final)
    psi = abs(psi)
    for i in range(len(psi)):
        psi[i] = psi[i]**2
    return psi

n = 2048
sigma = 3
matrix = np.zeros((n,n))

initial = np.zeros((n,1))
initial[n//2-1] = 1

J = 1
for i in range(n):
    for j in range(i,n,1):
        if i != j:
            u_ij = 2*(random.random()-0.5)
            #distance = min(abs(j-i),i+(n-j))
            #r_ij = (n/np.pi)*np.sin(np.pi*abs(distance)/n)
            r_ij = abs(j-i)
            t_ij = (J * u_ij)/((r_ij)**sigma)
            matrix[i,j] = t_ij
            matrix[j,i] = t_ij
t_array = [1,2,10]
for t in t_array:
    n_array = np.linspace(1,n,n)
    density = prob_density(wavefunc_t(matrix,n,initial,t))
    plt.plot(n_array,density)
    plt.xlim(800,1300)
plt.yscale('log')
plt.show()

