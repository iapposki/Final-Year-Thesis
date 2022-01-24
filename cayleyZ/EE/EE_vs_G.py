import numpy as np
import matplotlib.pyplot as plt
from cayley_gen import *

branch_indices = np.array([])
G = 5
Z = 3




#g = np.arange(2,11)     # number of generations
g = np.array([3,5,7,9])
Z = 3                   # coordination number
nf = 1                  # number of particles
EE_array = np.array([])
np.seterr(all='ignore')

for G in g:
    nA = int((nodes_till_generation(G,Z)-1)/Z)
    cayley = cayley_gen(Z,G)

    eigval_cayley,eigvec_cayley = np.linalg.eigh(cayley)

    c = np.array([])

    for i in range(G+1):
        if i == 2:
            branch_indices = np.array([2])
        elif i > 2:
            # print(i)
            for j in range(int(nodes_in_generation(i, Z) / Z)):
                branch_indices = np.append(branch_indices, [nodes_till_generation(i - 1, Z) + j + 1])
    print(G)
    nf = int(nodes_till_generation(G,Z)/2)
    #nf = 1

    for j in branch_indices:
        for i in branch_indices:
            add_this = 0
            for k in range(nf):
                add_this = add_this + eigvec_cayley[int(i), k] * (eigvec_cayley[int(j), k])

            #print(add_this)
            c = np.append(c, add_this)

    c = np.reshape(c,(len(branch_indices),len(branch_indices)))
    print(G)
    eigval,eigvec = np.linalg.eig(c)
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
plt.xscale('log')
plt.xlabel("Generations")
plt.ylabel("EE")
plt.title("Entanglement Entropy as a function of generations")
plt.plot(g,(EE_array))
plt.grid()
plt.savefig("test")
#np.savetxt("g.txt", g, delimiter=",")
np.savetxt("EE.txt", EE_array, delimiter=",")