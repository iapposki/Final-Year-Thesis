import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
import cmath


"""data_1 = np.loadtxt('EE_array_averaged_w=0.0.txt',dtype=float)
data_2 = np.loadtxt('EE_array_averaged_w=0.25.txt',dtype=float)
data_3 = np.loadtxt('EE_array_averaged_w=0.5.txt',dtype=float)
data_4 = np.loadtxt('EE_array_averaged_w=0.75.txt',dtype=float)
data_5 = np.loadtxt('EE_array_averaged_w=1.0.txt',dtype=float)
data_6 = np.loadtxt('EE_array_averaged_w=2.0.txt',dtype=float)

data_0 = np.loadtxt('n.txt',dtype=float)

data_1_err = np.loadtxt('err_w=0.0.txt',dtype=float)
data_2_err = np.loadtxt('err_w=0.25.txt',dtype=float)
data_3_err = np.loadtxt('err_w=0.5.txt',dtype=float)
data_4_err = np.loadtxt('err_w=0.75.txt',dtype=float)
data_5_err = np.loadtxt('err_w=1.0.txt',dtype=float)
data_6_err = np.loadtxt('err_w=2.0.txt',dtype=float)


plt.errorbar(data_0,data_1,yerr=data_1_err,label='0')
plt.errorbar(data_0,data_2,yerr=data_2_err,label='0.25')
plt.errorbar(data_0,data_3,yerr=data_3_err,label='0.50')
plt.errorbar(data_0,data_4,yerr=data_4_err,label='0.75')
plt.errorbar(data_0,data_5,yerr=data_5_err,label='1.0')
plt.errorbar(data_0,data_6,yerr=data_6_err,label='2.0')
plt.xscale('log')
plt.yticks([0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5])
plt.legend(title='W')
plt.xlabel("L")
plt.ylabel("EE")
plt.savefig('EE_vs_L.pdf')
plt.show()
"""

print(cmath.exp(1j*2))
