import numpy as np
import matplotlib.pyplot as plt
from cayley_gen import *

G = 5
Z = 3




#g = np.arange(2,11)     # number of generations
Gc = 9
g = np.array([9])
d = {2:128,3:64,4:32,5:16,6:8,7:4,8:2,9:1}
Z = 3                   # coordination number
nf = 1                  # number of particles
EE_array = np.array([])
np.seterr(all='ignore')
b_i = np.array([])



for G in g:
    #sub_sys = np.linspace(0, int((nodes_till_gen_bond_center(G) - 1)/2), int((nodes_till_gen_bond_center(G))/2))
    #sub_sys = np.array([63,127,128])
    cayley = cayley_gen_bond_center(Gc)
    print(len(cayley[0,:]))
    if G == 9:
        b = np.linspace(1, 255, 255,dtype='int')
    else:
        a = np.array([d[G]])
    #b = (neighbours_in_nearest_higher_generation_till_end(a, cayley) - 1)
    b = np.sort(b)
    sub_sys = b
    print(sub_sys)

    print("started")
    eigval_cayley,eigvec_cayley = np.linalg.eigh(cayley)
    print("stopped")
    c = np.array([])

    count = nodes_till_gen_bond_center(Gc)
    branch_indices = np.linspace(0,int((count-1)/2),int(count/2))
    #print(branch_indices)
    #branch_indices = neighbours_in_nearest_higher_generation_till_end(np.array([nodes_till_generation((Gc-G+1)-1,Z)+1]),cayley)
    b_i = np.append(b_i,len(sub_sys))
    print(G)
    nf = int(nodes_till_gen_bond_center(G)/2)
    #nf = 1

    for j in range(count):
        for i in range(count):
            add_this = 0
            for k in range(nf):
                add_this = add_this + eigvec_cayley[int(i), k] * (eigvec_cayley[int(j), k])

            #print(add_this)
            c = np.append(c, add_this)

    c = np.reshape(c,(count,count))
    #print(c)

    c_A = np.array([])
    for i in sub_sys:
        for j in sub_sys:
            c_A = np.append(c_A,c[int(i),int(j)])

    c_A = np.reshape(c_A,(len(sub_sys),len(sub_sys)))
    eigval,eigvec = np.linalg.eig(c_A)
    eigval = np.real(eigval)
    #print(eigval)

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
#plt.xscale('log')
plt.xlabel("Number of Lattices")
plt.ylabel("EE")
plt.title("Entanglement Entropy as a function of Number of Lattices")
plt.plot(b_i,(EE_array))
print(EE_array)
plt.grid()
plt.savefig("test")
np.savetxt('test.txt',EE_array,delimiter=',')

#np.savetxt("g.txt", g, delimiter=",")
#np.savetxt("EE.txt", EE_array, delimiter=",")

