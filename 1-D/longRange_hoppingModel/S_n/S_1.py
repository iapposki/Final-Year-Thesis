import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt
import math
import random
import scipy.linalg
import warnings
import cmath


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

def matrix_exp(power_matrix,t):
    #eigval,eigvec = scipy.linalg.eigh(power_matrix,overwrite_a=True,check_finite=True,turbo=True)
    eigval,eigvec = np.linalg.eig(power_matrix)
    #print(eigval)
    eigvec_inv = np.linalg.inv(eigvec)
    d = np.zeros((len(eigval),len(eigval)))
    d = d.astype("cdouble")
    for i in range(len(eigval)):
        d[i,i] = cmath.exp(1j*eigval[i]*t)
    res = np.matmul(eigvec,d)
    res = np.matmul(res,eigvec_inv)
    return res

def c_time(c_ini,H,t):
    a = matrix_exp(H,t)
    b = matrix_exp(H,-t)
    res = np.matmul(a,c_ini)
    res = np.matmul(res,b)
    #res = np.matmul(a,b)
    #print(a,b)
    return res

# L = number of lattice sites
#W = np.array([0])
t_array = np.logspace(-2,6,20)
L = 1024
nA = L//2
sigma_array = np.array([0.5,1.3,2.0,3.0])
#n =10**np.array([1,1.45,1.83,2,2.5,2.6,2.9])
#print(n)
numpy.seterr(all='ignore')
for sigma in sigma_array:
    EE_array_whole = np.array([])
    err = np.array([])
    print(f"doing for sigma = {sigma}")
    T = 20
    J = 1
    for x in range(T):
        EE_array = np.array([])
        print(x)
        for t in t_array:
            #nf = 1
            #nA = L     # Lattices in system A
            H = np.zeros((int(L), int(L)))
            initial = np.zeros((L,1))
            nF = 0
            for i in range(L):
                if i%2 == 0:
                    initial[i] = 1
                    nF += 1

            c_ini = np.zeros((L,L))
            for i in range(L):
                if i%2 != 0:
                    c_ini[i,i] = 1

            for i in range(L):
                for j in range(i,L,1):
                    if i != j:
                        u_ij = 2*(random.random()-0.5)
                        r_ij = abs(j-i)
                        t_ij = (J * u_ij)/((r_ij)**sigma)
                        H[i,j] = t_ij
                        H[j,i] = t_ij

            c_t = c_time(c_ini,H,t)
            c_t = c_t[:nF,:nF]

            eigval_c_t,eigvec_c_t = scipy.linalg.eigh(c_t,overwrite_a=False,check_finite=True,turbo=True)
            EE = 0
            for i in range(len(eigval_c_t)):
                p_A = eigval_c_t[i]
                if p_A == 1 or p_A == 0:
                    print("oof")
                else:
                    p_B = 1 - p_A
                    EE += - (p_A * np.log(abs(p_A))) - (p_B * np.log(abs(p_B)))

            EE_array = np.append(EE_array,np.real(EE))
        if x == 0:
            EE_array_whole = EE_array.copy()
        else:
            EE_array_whole = np.vstack([EE_array_whole, EE_array])

    for i in range(len(EE_array_whole[0, :])):
        err = np.append(err, np.std(EE_array_whole[:, i]))

    EE_array_whole_averaged = np.mean(EE_array_whole,axis=0)
    print(EE_array_whole_averaged)

    plt.errorbar(t_array, EE_array_whole_averaged, yerr=err, label=f"{sigma}")
    np.savetxt(f"EE_array_averaged_signa={sigma}_2.txt",EE_array_whole_averaged, delimiter=',')
    np.savetxt(f"err_sigma={sigma}_2.txt", err, delimiter=',')
#plt.plot(n,EE_array)
np.savetxt('t_array_2.txt', t_array, delimiter=',')

plt.legend(title="Sigma")
plt.xlabel("Time")
plt.ylabel("S_1")
plt.xscale("log")
plt.grid()
plt.savefig("test_2.pdf")
plt.show()
