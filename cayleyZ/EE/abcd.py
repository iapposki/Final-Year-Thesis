import numpy as np
import matplotlib.pyplot as plt
from cayley_gen import *

branch_indices = np.array([])
G = 5
Z = 3




#g = np.arange(2,11)     # number of generations
g = np.array([2,3,4,5,6])
d = {2:47,3:23,4:11,5:5,6:2}
Gc = 7
Z = 3                   # coordination number
nf = 1                  # number of particles
EE_array = np.array([])
np.seterr(all='ignore')
b_i = np.array([])

for G in g:
    nA = int((nodes_till_generation(G,Z)-1)/Z)
    cayley = cayley_gen(Z,Gc)
    np.savetxt('cayley.txt',cayley,delimiter=',',fmt='%d')
    print("started")
    eigval_cayley,eigvec_cayley = np.linalg.eigh(cayley)
    print("stopped")
    c = np.array([])

    """if G == 2:
        branch_indices = np.array([23,47,48])
    elif G == 3:
        branch_indices = np.array([11,23,24,47,48,49,50])
    elif G == 4:
        branch_indices = np.array([5,11,12,23,24,25,26,47,48,49,50,51,52,53,54])
    elif G == 5:
        branch_indices = np.array([2,5,6,11,12,13,14,23,24,25,26,27,28,29,30,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62])"""
    a = np.array([d[G]])
    b = neighbours_in_nearest_higher_generation_till_end(a, cayley) - 1
    b = np.sort(b)
    branch_indices = b
    b_i = np.append(b_i,len(branch_indices))
    print(G)
    nf = int(nodes_till_generation(Gc,Z)/2)
    #nf = 1

    for j in range(len(cayley[0,:])):
        for i in range(len(cayley[0,:])):
            add_this = 0
            for k in range(nf):
                add_this = add_this + eigvec_cayley[int(i), k] * (eigvec_cayley[int(j), k])

            #print(add_this)
            c = np.append(c, add_this)

    c = np.reshape(c,(len(cayley),len(cayley)))
    c_A = np.array([])
    for i in branch_indices:
        for j in branch_indices:
            c_A = np.append(c_A, c[int(i), int(j)])

    c_A = np.reshape(c_A, (len(branch_indices), len(branch_indices)))

    eigval,eigvec = np.linalg.eig(c_A)
    eigval = np.real(eigval)

    EE = 0
    #print(eigval)
    for i in range(len(eigval)):
        add_it = (abs(eigval[i])*np.log(abs(eigval[i])) + abs(1 - eigval[i])*np.log(abs(1-eigval[i])))
        #print(f"{add_it} and {i}")
        if add_it == np.inf:
            print("oof")
        elif add_it == -np.inf:
            print("foo")
        elif np.isnan(add_it) == True:
            #print("flip")
            True
        else:
            EE = EE - add_it
    EE_array = np.append(EE_array, np.real(EE))
print(b_i)
print(EE_array)
plt.plot(g,EE_array)
plt.savefig('test.png')
#plt.xscale('log')
#plt.xlabel("Generations")
#plt.ylabel("EE")
#plt.title("Entanglement Entropy as a function of generations")
#plt.plot(b_i,(EE_array))
#plt.grid()
#plt.savefig("test")
#np.savetxt("g.txt", g, delimiter=",")
#np.savetxt("EE.txt", EE_array, delimiter=",")