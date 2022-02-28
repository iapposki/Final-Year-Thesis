import numpy as np
import numpy.linalg
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


# L = number of lattice sites
#W = np.array([0])
L_size = np.array([4,8,16,32,64,128,512,600,700,800,900,1024])
L = 2048
sigma_array = np.array([0.5,1.3,2.0])
#n =10**np.array([1,1.45,1.83,2,2.5,2.6,2.9])
#print(n)
numpy.seterr(all='ignore')
for sigma in sigma_array:
    EE_array_whole = np.array([])
    err = np.array([])
    print(f"doing for sigma = {sigma}")
    T = 10
    J = 1
    for t in range(T):
        EE_array = np.array([])
        print(t)
        for nA in L_size:
            #nf = 1
            #nA = L     # Lattices in system A
            H = np.zeros((int(L), int(L)))
            for i in range(L):
                for j in range(i,L,1):
                    if i != j:
                        u_ij = 2*(random.random()-0.5)
                        r_ij = abs(j-i)
                        t_ij = (J * u_ij)/((r_ij)**sigma)
                        H[i,j] = t_ij
                        H[j,i] = t_ij

            eigval_H,eigvec_H = scipy.linalg.eigh(H,overwrite_a=False,check_finite=True,turbo=True)
            _EE = np.array([])
            for i in range(len(eigvec_H[0,:])):
                prob = eigvec_H[:,i]**2
                p_A = sum(prob[:nA])
                p_B = 1 - p_A
                EE = - (p_A * np.log(p_A)) - (p_B * np.log(p_B))
                _EE = np.append(_EE,EE)
            EE = sum(_EE)/len(_EE)
            EE_array = np.append(EE_array,np.real(EE))
        if t == 0:
            EE_array_whole = EE_array.copy()
        else:
            EE_array_whole = np.vstack([EE_array_whole, EE_array])

    for i in range(len(EE_array_whole[0, :])):
        err = np.append(err, np.std(EE_array_whole[:, i]))

    EE_array_whole_averaged = np.mean(EE_array_whole,axis=0)
    print(EE_array_whole_averaged)

    plt.errorbar(L_size, EE_array_whole_averaged, yerr=err, label=f"{sigma}")
    np.savetxt(f"EE_array_averaged_signa={sigma}.txt",EE_array_whole_averaged, delimiter=',')
    np.savetxt(f"err_sigma={sigma}.txt", err, delimiter=',')
#plt.plot(n,EE_array)
np.savetxt('L_size.txt', L_size, delimiter=',')

plt.legend(title="Sigma")
plt.xlabel("L")
plt.ylabel("EE")
plt.grid()
plt.savefig("test.pdf")
plt.show()
