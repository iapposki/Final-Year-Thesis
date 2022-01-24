import numpy as np
import matplotlib.pyplot as plt
import random
import warnings
import cmath
warnings.filterwarnings('ignore')


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
n = 710
t = 100
#sigma_array = np.linspace(0,16,9)
sigma_array = np.array([0.5,2,4,8])
print(sigma_array)
t_array = np.logspace(0,4,10)
t_array = t_array.tolist()
t_array = [0] + t_array

for sigma in sigma_array:
    avgr_array = np.array([])
    print(sigma , "------------------")
    count_arr = np.array([])
    for t in t_array:
        print(t)
        loop = 8
        count_avg = 0
        for l in range(loop):
            matrix = np.zeros((n,n))
            initial = np.zeros((n,1))
            nF = 0
            for i in range(n):
                if i%2 == 0:
                    initial[i] = 1
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
                        r_ij = abs(j-i)
                        t_ij = (J * u_ij)/((r_ij)**sigma)
                        matrix[i,j] = t_ij
                        matrix[j,i] = t_ij

            c_t = c_time(c_ini,matrix,t)
            c_t = c_t[:nF,:nF]
            sorted_eigval = np.linalg.eigvals(c_t)
            sorted_eigval = np.real(sorted_eigval)
            """for i in range(len(sorted_eigval)):
                #print(sorted_eigval[i],np.log((1-sorted_eigval[i]/sorted_eigval[i])))
                #temp = np.log(abs((1-sorted_eigval[i])/sorted_eigval[i]))
                temp = abs(sorted_eigval[i]/(1-sorted_eigval[i]))
                sorted_eigval[i] = temp"""
            count = 0
            for i in sorted_eigval:
                if (10**(-10))< i < (1 - 10**(-8)):
                    count += 1
            count_avg += count
        count_avg = count_avg/loop
        count_arr = np.append(count_arr,count_avg)
    plt.plot(t_array,count_arr,label=f"{sigma}")
plt.xlabel("time (t)")
plt.ylabel("S_0")
plt.legend(title='sigma', loc='upper left',bbox_to_anchor=(1, 0.75))
plt.xscale("log")
plt.tight_layout()
plt.savefig(f"s_0_temp(1).pdf",bbox_inches="tight")
plt.show()
