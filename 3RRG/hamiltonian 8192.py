import sys
import numpy as np
from sympy import *
import scipy.sparse.linalg
import random
import time
import matplotlib.pyplot as plt



def rand_cyclic(n):
    """
    :param n: Number of vertices of random cyclic graph
    :return: a : array, has order of vertices in random cyclic graph; b : Adjascency matrix of the random cyclic graph
    """

    # Making random cyclic graph containing all vertices. Array 'a' has the order of the vertices in the random cyclic graph
    a = np.array([], int)
    for k in range(n):
        a = np.append(a,k+1)
    a = np.random.choice(a,n,replace=False)
    #print(a)



    # Entries of matrix 'b' are 1/0, indicating if i-j vertex has an edge between them or not respectively.
    b = np.zeros((n,n))
    for k in range(n):
        i = a[k-1] - 1
        j = a[k] -1
        b[i,j] = b[i,j] + 1
        b[j,i] = b[j,i] + 1
    #print(b)

    return a,b

def three_rrg(n):
    """
    :param n: Number of vertices
    :return: The matrix b, the indices of which indicates if i-j vertex has an edge between them or not (1 or 0 respectively)
    """

    # We are making a random 3-regular graph

    a,b = rand_cyclic(n)

    m = n/2
    # Now we add edges such that the previously 2-regular graph is now a 3-regular graph.
    while len(a) > 0:
        l = len(a)
        i = a[0] - 1
        r = random.randint(1,l) - 1
        j = a[r] - 1

        while i == j or b[i,j]==1 :
            r = random.randint(1,l) - 1
            j = a[r] - 1
            m = m -1
            if m < 1 :
                break

        if m < 1:
            print("recursive ")
            three_rrg(n)
            a = np.array([])
        else :
            b[i,j] = b[i,j] + 1
            b[j,i] = b[j,i] + 1
            a = np.delete(a,0,axis=None)
            a = np.delete(a,r-1,axis=None)
        #print(b)            # The required random 3-regular graph
        #print(np.sum(b,axis=0))



    return b


begin = time.time()

def middle_values(test_list,K) :
    """
    :param a: array which needs to be sliced
    :param k: number of middle values which are needed
    :return: array which has k middle values of 'a'
    """

    # printing original list
    #print("The original list is : " + str(test_list))


    # computing strt, and end index
    strt_idx = (len(test_list) // 2) - (K // 2)
    end_idx = (len(test_list) // 2) + (K // 2)

    # slicing extracting middle elements
    res = test_list[strt_idx: end_idx + 1]

    return res


# Convert the Hamiltonian in matrix(sparse) and find viable method for finding eigenvalues.
# We already have hamiltonian for spinless particles hopping over a RRG with connectivity three in a potential disorder.

avgrrange = np.array([])           # Extent of disorder
n = 8192
avgavgrrange = np.array([])

for k in range(100):
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


        #eigval = np.linalg.eigvalsh(H)
        eigval = scipy.sparse.linalg.eigsh(H, k = 1024,sigma=0, return_eigenvectors=False,mode='normal')
        #print("")
        #print(eigvec)
        sorted_eigval = np.sort(eigval)
        #print(sorted_eigval)
        #sorted_eigval = middle_values(sorted_eigval,int(n/8))
        print(k)
        #print(sorted_eigval)
        level_spacing = np.array([])

        for i in range(1023):
            level_spacing = np.append(level_spacing,sorted_eigval[i+1]-sorted_eigval[i])
        #print(level_spacing)

        rn = np.array([])

        for i in range(1022) :
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
plt.title(' Mean adjacent gap ratio \'r\' as a function of disorder \'W\' at N = 8192')
plt.savefig("N=8192.png")
np.savetxt("N=8192.csv", finalr, delimiter=",")
np.savetxt("W_for_N=8192.csv", W, delimiter=",")
#plt.show()



end = time.time()
print("Time taken is " + str(end-begin) + " seconds")
