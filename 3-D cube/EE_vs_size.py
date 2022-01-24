import numpy as np
import scipy.sparse.linalg
import matplotlib.pyplot as plt



def cubic_adj(n):
    """
    :param n: number of lattice sites on the edge of cubic lattice
    :return: adjacency matrix of cubic lattice
    """

    adj = np.zeros((n**3,n**3))

    for i in range(n):
        for j in range(n):
            for k in range(n - 1):
                adj[n*n*i + n*j + k, n*n*i + n*j + k + 1] = 1
                adj[n*n*i + n*j + k + 1, n*n*i + n*j + k] = 1

    for i in range(n):
        for j in range(n - 1):
            for k in range(n):
                adj[n*n*i + n*j + k, n*n*i + n*j + k + n] = 1
                adj[n*n*i + n*j + k + n, n*n*i + n*j + k] = 1

    for i in range(n - 1):
        for j in range(n):
            for k in range(n):
                adj[n*n*i + n*j + k, n*n*i + n*j + k + n*n] = 1
                adj[n*n*i + n*j + k + n*n, n*n*i + n*j + k] = 1
    return adj


W = np.array([0,2,4,6,8,10,12,14,16,18,20])
n = np.array([3,4,5,6,7,8,9,10])
np.seterr(all='ignore')

for w in W:
    EE_array_whole = np.array([])
    err = np.array([])
    print(f"doing for w = {w}")
    if w == 0:
        T = 2
    else:
        T = 10
    for t in range(T):
        EE_array = np.array([])
        print(t)
        for L in n:

            print(L)
            nf = int(L**3 / 2)     # number of fermions in the lattice
            nA = int(L**3 / 2)     # Lattices in system A
            H = cubic_adj(int(L))

            for i in range(int(L**3)):
                H[i,i] = w*(np.random.random()-0.5)

            eigval_H,eigvec_H = np.linalg.eigh(H)

            c = np.zeros((nA,nA))

            for j in range(nA):
                for i in range(nA):
                    for k in range(nf):
                        c[i,j] = c[i,j] + eigvec_H[i,k]*(eigvec_H[j,k])

            #print(c)
            eigval,eigvec = np.linalg.eig(c)
            eigval = np.real(eigval)

            EE = 0
            #print(eigval)
            for i in range(len(eigval)):
                """if eigval[i] >= 1:
                    True
                elif eigval[i] <= 0:
                    True
                else:"""
                """print(np.log(abs(eigval[i])))
                print(np.log(abs(1-eigval[i])))"""

                add_it = (abs(eigval[i])*np.log(abs(eigval[i])) + abs(1 - eigval[i])*np.log(abs(1-eigval[i])))
                #print(f"{add_it} and {i}")
                if add_it == np.inf:
                    print("oof")
                elif add_it == -np.inf:
                    print("foo")
                elif np.isnan(add_it) == True:
                    #print("flip")
                    True
                else:
                    EE = EE - add_it

            EE_array = np.append(EE_array,np.real(EE))
        if t == 0:
            EE_array_whole = EE_array.copy()
        else:
            #print(len(EE_array))
            EE_array_whole = np.vstack([EE_array_whole, EE_array])
            #print(len(EE_array_whole[0,:]))

    for i in range(len(EE_array_whole[0, :])):
        err = np.append(err, np.std(EE_array_whole[:, i])/4)

    EE_array_whole_averaged = np.mean(EE_array_whole,axis=0)
    print(EE_array_whole_averaged)

    plt.errorbar(n**3, EE_array_whole_averaged, yerr=err, label=f"{w}")
    np.savetxt(f"EE_array_averaged_w={w}.txt",EE_array_whole_averaged, delimiter=',')
    np.savetxt(f"err_w={w}.txt", err, delimiter=',')
#plt.plot(n,EE_array)
np.savetxt('n.txt', n, delimiter=',')

plt.legend(title="W")
plt.xlabel("System Size")
plt.ylabel("EE")
plt.xscale('log')
#plt.grid()
plt.savefig("test")
plt.show()
