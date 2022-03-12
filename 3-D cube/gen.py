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

def rvW(W_min,W_max,delta,adj,loop,n):
    """
    :param W_min: Minimum strength of disorder
    :param W_max: Maximum strength of disorder
    :param delta: Step size of increasing disorder
    :param loop: Number of time r is averaged over
    :param adj: Adjacency matrix of the cubic lattice with edge length n
    :param n: edge length of cubic lattice
    :return: Array of W, r (Ensemble averaged ratio) and error (standard deviation scaled as 1/10 of original)
    """
    total_vert = n**3     #Total vertices in the cubic lattice
    avgrrange = np.array([])
    avgavgrrange = np.array([])
    err = np.array([])
    for k in range(loop):
        #print("ok - " + str(k))
        #print("ok1 - " + str(k))
        W = W_min - delta
        # print(np.linalg.eig(H))
        for j in range(int((W_max-W_min)/delta)):
            W = W + delta
            for i in range(total_vert) :
                w = W*(np.random.random() - 0.5)
                adj[i,i] = w

            l = n**3
            #eigval = scipy.sparse.linalg.eigs(adj, k = l, sigma=0, return_eigenvectors=False,maxiter=1000)
            eigval,eigvec = scipy.linalg.eigh(adj,overwrite_a=False,check_finite=True,turbo=True)
            sorted_eigval = np.real(np.sort(eigval))
            level_spacing = np.array([])
            for i in range(l-1):
                level_spacing = np.append(level_spacing,sorted_eigval[i+1]-sorted_eigval[i])
            #print(level_spacing)

            rn = np.array([])

            for i in range(l-2) :
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
        err = np.append(err,np.std(avgavgrrange[:,i])/4)
    finalr = np.mean(avgavgrrange,axis=0)
    W = np.linspace(W_min,W_max,int((W_max-W_min)/delta))
    return W,finalr,err

M = np.array([3,4,5,8])
for m in M:
    if m ==3:
        W,finalr,err = rvW(1,40,1,cubic_adj(m),1000,m)
    elif m == 4:
        W,finalr,err = rvW(1,40,1,cubic_adj(m),800,m)
    elif m == 5:
        W,finalr,err = rvW(1,40,1,cubic_adj(m),500,m)
    else:
        W,finalr,err = rvW(1,40,1,cubic_adj(m),100,m)
    print(m)
    plt.errorbar(W,finalr,yerr=err,label=f'{m}')
plt.xlabel('W')
plt.ylabel('r')
plt.ylim(0.28,0.6)
plt.legend(title='Edge Length')
plt.grid()
plt.savefig('r_vs_W.pdf')
plt.show()
