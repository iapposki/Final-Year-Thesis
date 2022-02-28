import numpy as np
import matplotlib.pyplot as plt
import random
import scipy.linalg

n = 1024
tau_range = np.logspace(-2,4,100)
sigma = np.array([0.1,0.3,0.5,1.3,2.0])
sigma_count = 0.5
matrix = np.zeros((n,n))
J = 1

for sigma_count in sigma:
    g = np.array([])
    err = np.array([])
    count = 0
    for tau in tau_range:
        T = 100
        g_temp = np.array([])
        for t in range(T):
            for i in range(n):
                for j in range(i+1,n,1):
                    u_ij = 2*(random.random()-0.5)
                    #distance = min(abs(j-i),i+(n-j))
                    #r_ij = (n/np.pi)*np.sin(np.pi*abs(distance)/n)
                    r_ij = float(abs(j-i))
                    t_ij = (J * u_ij)/((r_ij)**sigma_count)
                    matrix[i,j] = t_ij
                    matrix[j,i] = t_ij
            eigval = np.linalg.eigvalsh(matrix)
            summation = 0
            for i in range(len(eigval)):
                for j in range(len(eigval)):
                    #summation += np.cos(tau*(eigval[i]-eigval[j])) - j*np.sin(tau*(eigval[i]-eigval[j]))
                    summation += np.exp(-1j*tau*(eigval[i]-eigval[j]))
            summation = np.abs(summation)
            g_temp = np.append(g_temp,summation)
        g_std = np.std(g_temp)
        err = np.append(err,g_std)
        g_temp = sum(g_temp)/len(g_temp)
        g = np.append(g,g_temp)
        count += 1
        print(count)
    #plt.plot(tau_range,g,label=f"{sigma_count}")
    plt.errorbar(tau_range, g, yerr=err/10, label=f"{sigma_count}")
    np.savetxt(f"sigma_{sigma_count}.txt", g, delimiter=',')
    np.savetxt(f"err_{err}.txt", err, delimiter=',')
np.savetxt(f"tau_range.txt", tau_range, delimiter=',')
plt.xlabel("tau")
plt.ylabel("g(tau)")
plt.legend(title="sigma")
plt.xscale('log')
plt.yscale('log')
plt.savefig("res3_temp.pdf",dpi=1920)
plt.show()


