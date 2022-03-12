import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg
import random
# import seaborn



n_array = [512,1024,2048]
sigma_initial = 0
sigma = np.linspace(sigma_initial,3,15)
for n in n_array:
    print(f"((({n})))")
    matrix = np.zeros((n,n))
    J = 1
    avgr_array_array = np.array([])
    err = np.array([])
    """for i in range(n):
        w = W*(random.random() - 0.5)
        matrix[i,i] = w"""
    for W in range(50):
        avgr_array = np.array([])
        for sigma_count in sigma:
            for i in range(n):
                for j in range(i,n,1):
                    if i != j:
                        u_ij = 2*(random.random()-0.5)
                        #distance = min(abs(j-i),i+(n-j))
                        #r_ij = (n/np.pi)*np.sin(np.pi*abs(distance)/n)
                        r_ij = abs(j-i)
                        t_ij = (J * u_ij)/((r_ij)**sigma_count)
                        matrix[i,j] = t_ij
                        matrix[j,i] = t_ij

            eigval,eigvec = scipy.linalg.eigh(matrix,overwrite_a=False,check_finite=True,turbo=True)
            sorted_eigval = np.sort(eigval)
            level_spacing = np.array([])
            for i in range(n-1):
                level_spacing = np.append(level_spacing,sorted_eigval[i+1]-sorted_eigval[i])

            rn = np.array([])

            for i in range(n-2) :
                rn = np.append(rn,min(level_spacing[i],level_spacing[i+1])/max(level_spacing[i],level_spacing[i+1]))

            avgr = sum(rn)/len(rn)
            avgr_array = np.append(avgr_array,avgr)

        if W == 0:
            avgr_array_array = avgr_array
        else:
            avgr_array_array = np.vstack((avgr_array_array,avgr_array))
        print(W)

    for i in range(len(avgr_array_array[0,:])):
        err = np.append(err,np.std(avgr_array_array[:,i])/3)
    # for i in range(len(avgr_array_array[0,:])):
    #     if i != 2 and i != 6 and i != 12:
    #         err[i] = 0
    # print(err)
    avgr_array_array = np.mean(avgr_array_array, axis=0)
    # plt.plot(sigma,avgr_array_array,'-.',label=f"{n}")
    plt.errorbar(sigma,avgr_array_array,ls="-.",yerr=err, fmt='o', markersize=2, capsize=5,label=f"{n}")
    
plt.ylabel('$r$')
plt.xlabel('$\sigma$')
plt.legend(title="N")

x = np.linspace(0,3,10)
y = x-x+0.386
plt.plot(x,y,"k--")
x = np.linspace(0,3,10)
y = x-x + 0.529
plt.plot(x,y,"k--")

#np.savetxt("error.csv", err, delimiter=",")
#np.savetxt(f"avgr_N={n}.csv", n, delimiter=",")
#np.savetxt("sigma.csv", sigma, delimiter=",")
plt.savefig('r_vs_sigma_temp.pdf')
plt.show()
