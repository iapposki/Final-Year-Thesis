from cayley_gen import nodes_till_generation,cayley_gen,prob_density
import numpy as np
import math
import scipy.sparse.linalg
import random
import time
import matplotlib.pyplot as plt
import sympy as sp

begin = time.time()

total_gen = 5
z = 3
W = 5
time_val = 501

total_vert = nodes_till_generation(total_gen,z)
initial = np.zeros((total_vert,1))
initial[0] = 1
############## Valid only for this initial, further modifications needed for a general initial in hist_weight !!!!!!!!!!
cayley = cayley_gen(z,total_gen)
def wavefunc_t(cayley,total_vert,initial,t):
        eigval,eigvec = scipy.linalg.eigh(cayley,overwrite_a=True,check_finite=True,turbo=True)
        #eigval,eigvec = np.linalg.eig(H)
        #eigval,eigvec = scipy.sparse.linalg.eigsh(H, k = int(n/2), sigma=0, return_eigenvectors=True)
        psi = np.zeros(total_vert)
        for i in range(total_vert):
            psi = psi + math.e**(-1j*t*eigval[i])*(np.dot(np.transpose(eigvec[:,i]),initial))*eigvec[:,i]
            #psi = psi + sp.exp((-1j*t*eigval[i]))*((sp.transpose(sp.Matrix(eigvec[:,i])).dot(sp.Matrix(initial))))*sp.Matrix(eigvec[:,i])
        return psi

def hist_weight(psi,total_gen,z):
        """
        :param psi: wavefunction amplitudes.
        :param total_gen: total generations.
        :param z: coordination number.
        :return: probability density apmlitudes as function of generations.
        """
        psi = prob_density(psi)
        hist = np.zeros((total_gen,1))

        for i in range(total_gen):
            if i == 0:
                hist[0] = i*psi[0]
            else:
                hist[i] = i*sum(psi[nodes_till_generation(i-1,z)+1:nodes_till_generation(i,z)+1])
        hist_weight_sum = sum(hist)
        return hist_weight_sum

def weightFuncOf_t_vs_Time(total_vert,total_gen,z,cayley,initial,total_time,W):
    # weight is the summation of probability magnitude at the lattice site times the generation of that lattice, center
    # lattice is given zero weight, i.e. more the spread, more the weight.
    for i in range(total_vert) :
        w = W*(np.random.random() - 0.5)
        cayley[i,i] = w
    cayley_cpy = cayley.copy()
    weight_func_t = np.array([])
    for i in range(total_time):
        weight_func_t = np.append(weight_func_t,hist_weight(wavefunc_t(cayley,total_vert,initial,i/10),total_gen,z))
    time_points = np.linspace(0,total_time-1,total_time)
    return weight_func_t,time_points




###############
# weight as a function of time
###############
"""weight_func_t,time_points = weightFuncOf_t_vs_Time(total_vert,total_gen,z,cayley,initial,time_val,W)
plt.plot(time_points,weight_func_t,'r-')
plt.grid()
plt.xlabel('Time')
plt.ylabel('Weight')
plt.title('Weight as a function of time')
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
print(f"mean of the weight function for W = {W} is {np.mean(weight_func_t)}, with standard deviation {np.std(weight_func_t)}")
plt.show()
"""


###############
# weight function average, standard deviation, standard deviation of standard deviation for different coordination
# number of cayley trees at long enough time.
###############
"""weightfunction_avg = np.array([])
std_avg = np.array([])
std_std_avg = np.array([])
for j in [3,4,5,6,7,8]:
    z = j
    temp1 = np.array([])
    temp2 = np.array([])
    print(z)
    for i in range(10):
        weight_func_t,time_points = weightFuncOf_t_vs_Time(total_vert,total_gen,z,cayley,initial,time_val,W)
        #plt.plot(time_points,weight_func_t,'r-')
        #plt.grid()
        #figManager = plt.get_current_fig_manager()
        #figManager.window.showMaximized()
        #print(f"mean of the weight function for W = {W} is {np.mean(weight_func_t)}, with standard deviation {np.std(weight_func_t)}")
        temp1 = np.append(temp1,np.mean(weight_func_t))
        temp2 = np.append(temp2,(np.std(weight_func_t))**2)
        #plt.show()
    weightfunction_avg = np.append(weightfunction_avg,np.mean(temp1))
    std_avg = np.append(std_avg,np.sqrt(np.mean(temp2)))
print(weightfunction_avg)
print(std_avg)
"""




##########################
# weight function of increasing disorder at long enough time
##########################

"""weightfunction_avg = np.array([])
std_avg = np.array([])
for j in range(51):
    W = j
    temp1 = np.array([])
    temp2 = np.array([])
    print(W)
    for i in range(20):
        weight_func_t,time_points = weightFuncOf_t_vs_Time(total_vert,total_gen,z,cayley,initial,time_val,W)
        #plt.plot(time_points,weight_func_t,'r-')
        #plt.grid()
        #figManager = plt.get_current_fig_manager()
        #figManager.window.showMaximized()
        #print(f"mean of the weight function for W = {W} is {np.mean(weight_func_t)}, with standard deviation {np.std(weight_func_t)}")
        temp1 = np.append(temp1,np.mean(weight_func_t))
        temp2 = np.append(temp2,np.std(weight_func_t))
        #plt.show()
    weightfunction_avg = np.append(weightfunction_avg,np.mean(temp1))
    std_avg = np.append(std_avg,np.mean(temp2))

W_range = np.linspace(0,50,51)

plt.errorbar(W_range,weightfunction_avg,yerr=std_avg)
plt.xlabel('W')
plt.ylabel('Averaged Weight Function')
plt.title('Averaged Weight Function as a function of W; Z = 3')
plt.grid()
np.savetxt("W_range.csv", W_range, delimiter=",")
np.savetxt("avg_std_dev_vs_W_z=3.csv", std_avg, delimiter=",")
plt.savefig('weight_func_vs_W_z=3.png', format='png', dpi=1920)
plt.show()"""



############################
# Multiple subplots for given
# disorder over a range
############################
time_val = 101
z = 3
W = 0
total_gen = 3
total_vert = nodes_till_generation(total_gen,z)
initial = np.zeros((total_vert,1))
initial[0] = 1
cayley = cayley_gen(z,total_gen)
plt.rcParams['axes.grid'] = True
for k in range(26):
    err = np.array([])
    for j in range(30):
        print(j)
        #plt.title(f'Weight plots; W = {W}; Z = {z}')
        weight_func_t,time_points = weightFuncOf_t_vs_Time(total_vert,total_gen,z,cayley,initial,time_val,W)
        if j == 0:
            weight_func_t_whole = weight_func_t.copy()
        else:
            weight_func_t_whole = np.vstack([weight_func_t_whole,weight_func_t])
    for i in range(len(weight_func_t_whole[0, :])):
        err = np.append(err, np.std(weight_func_t_whole[:, i]))
    weight_func_t_averaged = np.mean(weight_func_t_whole,axis=0)
    plt.errorbar(time_points/10, weight_func_t_averaged, yerr=err, color='red')
    plt.xlabel("Time")
    plt.ylabel("Weight")
    plt.ylim(0,3.2)
    #plt.tight_layout()
    plt.title(f"Averaged Weight vs Time; Z={z} W={W} in {total_gen} gens")
    plt.savefig(f"weight_func_vs_t_Z={z}_W={W}gen={total_gen}.png", format="png", dpi=1280)
    #plt.show()
    plt.clf()
    W = W + 2