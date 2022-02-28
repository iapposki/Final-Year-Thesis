import numpy as np
import matplotlib.pyplot as plt


data_1 = np.loadtxt('EE_array_averaged_w=0.txt',dtype=float)
data_2 = np.loadtxt('EE_array_averaged_w=2.txt',dtype=float)
data_3 = np.loadtxt('EE_array_averaged_w=4.txt',dtype=float)
data_4 = np.loadtxt('EE_array_averaged_w=6.txt',dtype=float)
data_5 = np.loadtxt('EE_array_averaged_w=8.txt',dtype=float)
data_6 = np.loadtxt('EE_array_averaged_w=10.txt',dtype=float)
data_7 = np.loadtxt('EE_array_averaged_w=12.txt',dtype=float)
data_8 = np.loadtxt('EE_array_averaged_w=14.txt',dtype=float)
data_9 = np.loadtxt('EE_array_averaged_w=16.txt',dtype=float)
data_10 = np.loadtxt('EE_array_averaged_w=18.txt',dtype=float)
data_11 = np.loadtxt('EE_array_averaged_w=20.txt',dtype=float)

data_0 = np.loadtxt('n.txt',dtype=float)
data_0 = data_0**2

data_1_err = np.loadtxt('err_w=0.txt',dtype=float)
data_2_err = np.loadtxt('err_w=2.txt',dtype=float)
data_3_err = np.loadtxt('err_w=4.txt',dtype=float)
data_4_err = np.loadtxt('err_w=6.txt',dtype=float)
data_5_err = np.loadtxt('err_w=8.txt',dtype=float)
data_6_err = np.loadtxt('err_w=10.txt',dtype=float)
data_7_err = np.loadtxt('err_w=12.txt',dtype=float)
data_8_err = np.loadtxt('err_w=14.txt',dtype=float)
data_9_err = np.loadtxt('err_w=16.txt',dtype=float)
data_10_err = np.loadtxt('err_w=18.txt',dtype=float)
data_11_err = np.loadtxt('err_w=20.txt',dtype=float)


plt.errorbar(data_0,data_1,yerr=data_1_err,label='0')
plt.errorbar(data_0,data_2,yerr=data_2_err,label='2')
plt.errorbar(data_0,data_3,yerr=data_3_err,label='4')
plt.errorbar(data_0,data_4,yerr=data_4_err,label='6')
plt.errorbar(data_0,data_5,yerr=data_5_err,label='8')
plt.errorbar(data_0,data_6,yerr=data_6_err,label='10')
plt.errorbar(data_0,data_7,yerr=data_7_err,label='12')
plt.errorbar(data_0,data_8,yerr=data_8_err,label='14')
plt.errorbar(data_0,data_9,yerr=data_9_err,label='16')
plt.errorbar(data_0,data_10,yerr=data_10_err,label='18')
plt.errorbar(data_0,data_11,yerr=data_11_err,label='20')

#plt.yticks([0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5])
plt.legend(title='W')
plt.xlabel("L")
plt.ylabel("EE")
plt.savefig('EE_vs_L_test.pdf')
plt.show()
