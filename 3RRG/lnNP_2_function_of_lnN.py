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

m = np.array([10000,5000,1000,200,100,50,20,10])
n = np.array([64,128,256,512,1024,2048,4096,8192])
W = np.array([0,1,2,3,5,9,10,12,14,18,25])
lnpn = np.array([])
for l in range(len(W)):
    for k in range(len(n)):
        time_val = 10
        final_all = np.zeros((n[k],1))

        initial = np.zeros((n[k],1))
        initial[int(n[k]/2)] = 1
        total_for_average = 0
        for j in range(int(m[k])):
            H = three_rrg(n[k])
            for i in range(n[k]) :
                w = W[l]*(random.random() - 0.5)
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

            final = wavefunc_t(H_cpy,n[k],initial,time_val)

            #final = final.evalf()
            final = abs(final)
            #print(final)

            for i in range(n[k]):
                final[i] = final[i]**2

            total = 0
            for i in range(n[k]):
                total = total + final[i]**2

            total_for_average = total_for_average + total*n[k]

        average = total_for_average/m[k]
        lnpn_val = math.log(average)
        lnpn = np.append(lnpn,[lnpn_val])
        #print(abs(Dagger(final).dot(final)).evalf())
        #print(final[0])

        #print(sp.transpose(sp.Matrix(final)).dot(sp.Matrix(final)))

    if l==0 :
        wholematrix = lnpn.copy()
    else :
        wholematrix = np.vstack([wholematrix,lnpn])
    lnpn = np.array([])
    checkpoint = time.time()
    print(f"{int((l+1)*100/len(W))}% completed; Time passed : {round(int(checkpoint-begin)/60.0,1)} minutes")

print(wholematrix)
for j in range(len(W)):
    plt.plot(np.log(n),wholematrix[j,:],label=f"{W[j]}")

np.savetxt("Datapoints_for_lnNPvslnN.csv", wholematrix, delimiter=",")
np.savetxt("lnN_for_lnNPvslnN.csv", np.log(n), delimiter=",")

plt.xlabel('ln N')
plt.ylabel('ln(N P)')
plt.legend(title='W',loc=2, bbox_to_anchor=(1.01, 1))
plt.grid()
plt.savefig('lnNP_vs_lnN.png', format='png', dpi=1920, bbox_inches="tight")
plt.show()

end = time.time()
print("Time taken is " + str(end-begin) + " seconds")
