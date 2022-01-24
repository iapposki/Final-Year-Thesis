import matplotlib.pyplot as plt
import numpy as np
from cayley_gen import *
import random
import scipy.linalg


cayley = cayley_gen(3,2)
eigval,eigvec = scipy.linalg.eigh(cayley,overwrite_a=True,check_finite=True,turbo=True)
print(eigvec)
