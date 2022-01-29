import numpy as np
import matplotlib.pyplot as plt
import cmath
import random
import scipy.linalg
import sympy as sp
import warnings
#warnings.filterwarnings('ignore')



def wavefunc_t(matrix,n,initial,t):
    #eigval,eigvec = scipy.linalg.eigh(matrix,overwrite_a=True,check_finite=True,turbo=True)
    eigval,eigvec = np.linalg.eig(matrix)
    #eigval,eigvec = scipy.sparse.linalg.eigsh(H, k = int(n/2), sigma=0, return_eigenvectors=True)
    psi = np.zeros(n)
    for i in range(n):
        psi = psi + cmath.e**(-1j*t*eigval[i])*(np.dot(np.transpose(eigvec[:,i]),initial))*eigvec[:,i]
        #psi = psi + sp.exp((-1j*t*eigval[i]))*((sp.transpose(sp.Matrix(eigvec[:,i])).dot(sp.Matrix(initial))))*sp.Matrix(eigvec[:,i])
    return psi

def matrix_exp(power_matrix,t):
    #eigval,eigvec = scipy.linalg.eigh(power_matrix,overwrite_a=True,check_finite=True,turbo=True)
    eigval,eigvec = np.linalg.eig(power_matrix)
    #print(eigval)
    eigvec_inv = np.transpose(eigvec)
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

n = 1024
t = 100
#sigma_array = np.linspace(0,16,9)
sigma_array = np.array([0.5,2.0,2.7,3.0,3.5,4.0,8.0])
print(sigma_array)
t_array = np.logspace(0,6,20)
for sigma in sigma_array:
    avgr_array = np.array([])
    print(sigma)
    for t in t_array:
        print(t)
        loop = 100
        avg_avgr = 0
        for l in range(loop):
            matrix = np.zeros((n,n))
            #initial = np.zeros((n,1))
            nF = 0
            for i in range(n):
                if i%2 == 0:
                    #initial[i] = 1
                    nF += 1

            c_ini = np.zeros((n,n))
            for i in range(n):
                if i%2 != 0:
                    c_ini[i,i] = 1

            J = 1
            for i in range(n):
                for j in range(i,n,1):
                    if i != j:
                        u_ij = 2*(random.random()-0.5)
                        r_ij = float(abs(j-i))

                        if pow(r_ij,sigma) == 0:
                            t_ij = 0
                            print("oof")
                        else:
                            t_ij = (J * u_ij)/(pow(r_ij,sigma))
                            #print(t_ij)
                        matrix[i,j] = t_ij
                        matrix[j,i] = t_ij
            c_t = c_time(c_ini,matrix,t)
            c_t = c_t[:nF,:nF]
            sorted_eigval = np.linalg.eigvals(c_t)
            #sorted_eigval = np.real(sorted_eigval)
            indices = []
            for i in range(len(sorted_eigval)):
                #print(sorted_eigval[i],np.log((1-sorted_eigval[i]/sorted_eigval[i])))
                #sorted_eigval[i] = float(sorted_eigval[i])
                if (sorted_eigval[i]) != 0 and (sorted_eigval[i]) != 1:
                    #print(sorted_eigval[i])
                    temp = (np.log(abs((1-(sorted_eigval[i]))/(sorted_eigval[i]))))
                    sorted_eigval[i] = float(temp)
                else:
                    indices.append(i)
            sorted_eigval = np.delete(sorted_eigval,indices)
                
            level_spacing = np.array([])
            sorted_eigval = np.sort(sorted_eigval)
            for i in range(len(sorted_eigval)-1):
                level_spacing = np.append(level_spacing,sorted_eigval[i+1]-sorted_eigval[i])
                #print(level_spacing[-1])

            rn = np.array([])

            for i in range(len(sorted_eigval)-2) :
                temp = min(level_spacing[i],level_spacing[i+1])/max(level_spacing[i],level_spacing[i+1])
                if not np.isnan(temp):
                    rn = np.append(rn,temp)

            avgr = sum(rn)/len(rn)
            avg_avgr += avgr
        avgr = avg_avgr/loop
        avgr_array = np.append(avgr_array,avgr)
    plt.plot(t_array,avgr_array,label=f"{sigma}")
    np.savetxt(f"sigma_{sigma}.txt", avgr_array, delimiter=',')
np.savetxt('t_array.txt', t_array, delimiter=',')
plt.xlabel("time (t)")
plt.ylabel("<r>")
plt.legend(title="sigma", loc='upper left',bbox_to_anchor=(1, 0.75))
plt.xscale("log")
plt.savefig("res.png", bbox_inches="tight")
#plt.show()
