from cayley_gen import *
import numpy as np

"""z = 3
total_gen = 4
cayley = cayley_gen(z,total_gen)

np.savetxt('test.txt', cayley,fmt='%d', delimiter=",")
"""

a = np.array([1,2,3,4,5])
b = np.array([-1,0])
print(min(a)/max(b))
