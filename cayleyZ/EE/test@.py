import numpy as np

from cayley_gen import *

"""cayley = cayley_gen(3,3)
print(len(cayley[0,:]))
print(nodes_till_generation(11,3))
#np.savetxt('cayley.txt',cayley,delimiter=',',fmt='%d')
a1 = np.array([63])
a2 = np.array([2])
#b = np.append(neighbours_in_nearest_higher_generation_till_end(a2,cayley) - 1,neighbours_in_nearest_higher_generation_till_end(a1,cayley) - 1)
b = neighbours_in_nearest_higher_generation_till_end(a2,cayley)
b = np.sort(b)
print(b)"""
profit = 0
nums = [8,7,6,5,4,5,6,7,8,9]

for i in range(len(nums) - 1):
    print(nums[i])
    print(nums[i+1:])
    print(nums[i] in nums[i + 1:])
    if nums[i] in nums[i + 1:]:
        print(True)
    else:
        print(False)