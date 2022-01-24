import numpy as np
import random

def rand_cyclic(n):
    """
    :param n: Number of vertices of random cyclic graph
    :return: a : array, has order of vertices in random cyclic graph; b : Adjascency matrix of the random cyclic graph
    """
    # Making random cyclic graph containing all vertices. Array 'a' has the order of the vertices in the random cyclic graph
    a = np.array([], int)
    for k in range(n):
        a = np.append(a,k+1)
    a = np.random.choice(a,n,replace=False)
    #print(a)



    # Entries of matrix 'b' are 1/0, indicating if i-j vertex has an edge between them or not respectively.
    b = np.zeros((n,n))
    for k in range(n):
        i = a[k-1] - 1
        j = a[k] -1
        b[i,j] = b[i,j] + 1
        b[j,i] = b[j,i] + 1
    #print(b)

    return a,b

def three_rrg(n):
    """
    :param n: Number of vertices
    :return: The matrix b, the indices of which indicates if i-j vertex has an edge between them or not (1 or 0 respectively)
    """
    # We are making a random 3-regular graph
    a,b = rand_cyclic(n)

    m = n/2
    # Now we add edges such that the previously 2-regular graph is now a 3-regular graph.
    while len(a) > 0:
        l = len(a)
        i = a[0] - 1
        r = random.randint(1,l) - 1
        j = a[r] - 1

        while i == j or b[i,j]==1 :
            r = random.randint(1,l) - 1
            j = a[r] - 1
            m = m -1
            if m < 1 :
                break

        if m < 1:
            three_rrg(n)
            a = np.array([])
        else :
            b[i,j] = b[i,j] + 1
            b[j,i] = b[j,i] + 1
            a = np.delete(a,0,axis=None)
            a = np.delete(a,r-1,axis=None)
        #print(b)            # The required random 3-regular graph
        #print(np.sum(b,axis=0))



    return b


def middle_values(test_list,K) :
    """
    :param a: array which needs to be sliced
    :param k: number of middle values which are needed
    :return: array which has k middle values of 'a'
    """

    # printing original list
    #print("The original list is : " + str(test_list))


    # computing strt, and end index
    strt_idx = (len(test_list) // 2) - (K // 2)
    end_idx = (len(test_list) // 2) + (K // 2)

    # slicing extracting middle elements
    res = test_list[strt_idx: end_idx + 1]

    return res
