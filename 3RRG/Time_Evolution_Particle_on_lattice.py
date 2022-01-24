import math
import sys
import numpy as np
from sympy import *
import scipy.sparse.linalg
from Making_graph_matrix import *
import random
import time
import matplotlib.pyplot as plt
import sympy as sp
from sympy import *
from sympy.physics.quantum import *
import scipy.linalg

begin = time.time()

n = 512
H = three_rrg(n)
W = 12
time_val = 100



final_all = np.zeros((n,1))
initial = np.zeros((n,1))
initial[int(n/2)] = 1


for i in range(n) :
    w = W*(random.random() - 0.5)
    H[i,i] = w
H_cpy = H.copy()
def wavefunc_t(H,n,initial,t):
    eigval,eigvec = scipy.linalg.eigh(H,overwrite_a=True,check_finite=True,turbo=True)
    #eigval,eigvec = np.linalg.eig(H)
    #eigval,eigvec = scipy.sparse.linalg.eigsh(H, k = int(n/2), sigma=0, return_eigenvectors=True)
    psi = np.zeros(n)
    for i in range(n):
        psi = psi + math.e**(-1j*t*eigval[i])*(np.dot(np.transpose(eigvec[:,i]),initial))*eigvec[:,i]
        #psi = psi + sp.exp((-1j*t*eigval[i]))*((sp.transpose(sp.Matrix(eigvec[:,i])).dot(sp.Matrix(initial))))*sp.Matrix(eigvec[:,i])
    return psi

final = wavefunc_t(H_cpy,n,initial,time_val)

#final = final.evalf()
final = abs(final)
#print(final)
for i in range(n):
    final[i] = final[i]**2

total = 0
for i in range(n):
    total = total + final[i]**2
print(n*total)

shape(final)
#print(abs(Dagger(final).dot(final)).evalf())
#print(final[0])

#print(sp.transpose(sp.Matrix(final)).dot(sp.Matrix(final)))

end = time.time()
print("Time taken is " + str(end-begin) + " seconds")


n = np.linspace(1,n,n)
plt.plot(n,final,'r-')
plt.show()



