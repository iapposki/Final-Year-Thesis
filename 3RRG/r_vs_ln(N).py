import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline


N64 = np.genfromtxt('N=64.csv', delimiter="\t")
N128 = np.genfromtxt('N=128.csv', delimiter="\t")
N256 = np.genfromtxt('N=256.csv', delimiter="\t")
N512 = np.genfromtxt('N=512.csv', delimiter="\t")
N1024 = np.genfromtxt('N=1024.csv', delimiter="\t")
N2048 = np.genfromtxt('N=2048.csv', delimiter="\t")
N4096 = np.genfromtxt('N=4096.csv', delimiter="\t")
N8192 = np.genfromtxt('N=8192.csv', delimiter="\t")
W = np.genfromtxt('W_for_N=64.csv', delimiter="\t")


logN = np.array([np.log(64),np.log(128),np.log(256),np.log(512),np.log(1024),np.log(2048),np.log(4096),np.log(8192)])

#print(np.transpose(N64))

r = np.vstack([N64,N128,N256,N512,N1024,N2048,N4096,N8192])
r = np.transpose(r)



plt.plot(logN,r[0,:],label='5',color='#0000ff')
plt.plot(logN,r[4,:],label='9',color='#0025dc')
plt.plot(logN,r[6,:],label='11',color='#0049b7')
plt.plot(logN,r[7,:],label='12',color='#006c93')
plt.plot(logN,r[8,:],label='13',color='#00906c')
plt.plot(logN,r[9,:],label='14',color='#00b448')
plt.plot(logN,r[12,:],label='17',color='#00da24')
plt.plot(logN,r[19,:],label='25',color='#00ff00')
plt.xlabel('ln N')
plt.ylabel('r')
plt.legend(title='W', loc='upper left',bbox_to_anchor=(1, 0.75))
plt.savefig("r_vs_lnN.png", format='png', dpi=1920, bbox_inches="tight")
plt.show()
