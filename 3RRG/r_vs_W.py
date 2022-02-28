import numpy as np
import matplotlib.pyplot as plt



N64 = np.genfromtxt('N=64.csv', delimiter="\t")
N128 = np.genfromtxt('N=128.csv', delimiter="\t")
N256 = np.genfromtxt('N=256.csv', delimiter="\t")
N512 = np.genfromtxt('N=512.csv', delimiter="\t")
N1024 = np.genfromtxt('N=1024.csv', delimiter="\t")
N2048 = np.genfromtxt('N=2048.csv', delimiter="\t")
N4096 = np.genfromtxt('N=4096.csv', delimiter="\t")
N8192 = np.genfromtxt('N=8192.csv', delimiter="\t")
W = np.genfromtxt('W_for_N=64.csv', delimiter="\t")



plt.plot(W,N64,"",W,N128,"",W,N256,"",W,N512,"",W,N1024,"",W,N2048,"",W,N4096,"",W,N8192,"")
plt.plot(W,N64,color='#0000ff',label='64')
plt.plot(W,N128,color='#0025dc',label='128')
plt.plot(W,N256,color='#0049b7',label='256')
plt.plot(W,N512,color='#006c93',label='512')
plt.plot(W,N1024,color='#00906c',label='1024')
plt.plot(W,N2048,color='#00b448',label='2048')
plt.plot(W,N4096,color='#00da24',label='4096')
plt.plot(W,N8192,color='#00ff00',label='8192')
plt.grid()
#plt.legend(["64","128","256","512","1024","2048","4096","8192"], title='N', loc ="best")
plt.legend(title='N',loc='best')
plt.xlabel('W')
plt.ylabel('r')
plt.savefig("r_vs_W.png", format='png', dpi=1920)
plt.show()
