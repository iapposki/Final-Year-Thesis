import sys
import numpy as np
from sympy import *
import scipy.sparse.linalg
from Making_graph_matrix import *
import random
import time
import matplotlib.pyplot as plt

begin = time.time()



# Convert the Hamiltonian in matrix(sparse) and find viable method for finding eigenvalues.
# We already have hamiltonian for spinless particles hopping over a RRG with connectivity three in a potential disorder.

avgrrange = np.array([])           # Extent of disorder
n = 256
avgavgrrange = np.array([])

for k in range(15000):
    print("ok - " + str(k))
    H = three_rrg(n)
    #print("ok1 - " + str(k))
    W = 5
    delta = 1
    # print(np.linalg.eig(H))
    for j in range(int(20/delta)):
        W = W + delta
        for i in range(n) :
            w = W*(random.random() - 0.5)
            H[i,i] = w
        #print(H)


        eigval = scipy.sparse.linalg.eigs(H, k = int(n/8), sigma=0, return_eigenvectors=False)
        eigval = np.real(eigval)
        #print(eigval)
        #print("")
        #print(eigvec)
        sorted_eigval = np.sort(eigval)
        #print(sorted_eigval)
        level_spacing = np.array([])
        for i in range(int(n/8)-1):
            level_spacing = np.append(level_spacing,sorted_eigval[i+1]-sorted_eigval[i])
        #print(level_spacing)

        rn = np.array([])

        for i in range(int(n/8)-2) :
            rn = np.append(rn,min(level_spacing[i],level_spacing[i+1])/max(level_spacing[i],level_spacing[i+1]))

        avgr = sum(rn)/len(rn)
        #print(avgr)
        avgrrange = np.append(avgrrange,avgr)
    if k == 0 :
        avgavgrrange = avgrrange.copy()
    else :
        avgavgrrange = np.vstack([avgavgrrange,avgrrange])
    avgrrange = np.array([])
finalr = np.mean(avgavgrrange,axis=0)





W = np.linspace(5,25,int(20/delta))



plt.plot(W,finalr,'r-')
plt.xlabel('W')
plt.ylabel('r')
plt.title(' Mean adjacent gap ratio \'r\' as a function of disorder \'W\' at N = 256')
plt.savefig("N=256.png")
np.savetxt("N=256.csv", finalr, delimiter=",")
np.savetxt("W_for_N=256.csv", W, delimiter=",")
plt.show()



end = time.time()
print("Time taken is " + str(end-begin) + " seconds")
