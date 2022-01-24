# this functions rounds the given number to needed significant numbers
def round_sig(x, sig) :
    from math import log10, floor
    return round(x, sig-int(floor(log10(abs(x))))-1)
def bisection_root(f, a, b, err) :
    an = a
    bn = b
    c = (a + b) / 2
    co = a
    cn = c
    while abs((cn - co) / cn) > err:
        co = cn
        if f(an,g) * f(cn) < 0:
            bn = cn
            cn = (bn + an) / 2
        elif f(cn) * f(bn) < 0:
            an = cn
            cn = (an + bn) / 2
        elif f(an) * f(bn) == 0:
            if f(an) == 0:
                print(str(an) + "is the root of the equation in given interval")
            else:
                print(str(bn) + "is the root of the equation in given interval")

    print(str(cn) + " is the root of the equation in the given interval for relative error = " + str(err))
def newton_rapson(g,x0,err) :
    """
    :param g: explicitfunction of 'x'(NOTE : make sure to use sympy for trigonometric and other functions)
    :param x0: initial guess
    :param err: required error
    :return: root corresponding to specifeid error
    """
    import math
    import sympy
    from sympy import diff, symbols
    import numpy
    xo = x0
    x = symbols('x')
    xn = x0 - float(g(x0) / (diff(g(x), x).subs(x, x0)))
    j = 0

    while abs((xn - xo) / xn) > err:
        xo = xn
        x = symbols('x')
        xn = xo - float(g(xo) / (diff(g(x), x).subs(x, xo)))
        if xn == 0 :
            break
        j += 1
    return xn
def gauss_el_roots(d) :
    """
    :param d:   'N x N+1' matrix
    :return:    'N x 1' matrix containing roots corresponding to matrix d
    """
    from numpy import zeros
    b = d.copy()
    n = len(b[:,0])
    c = zeros((n,1))
    i = 0
    while i < n :           # this part converts the matrix to row echelon form
        j = 0
        while j + i + 1 < n :
            if b[i, i] != 0 :
                b[j + i + 1, :] = b[j + i + 1, :] - b[i, :]*float(b[j + i + 1, i]/b[i, i])
                j += 1
        i += 1
    #print(b)
    j = 0
    while j < n :           # this part makes the leading entries of each row 1
        b[j,:] = b[j, :]/b[j, j]
        j += 1
    #print(b)
    j = 0
    i = 0
    while n - j - 1 > 0 :   # this part convert the matrix to reduced row echelon form
        i = 0
        while n - i - j - 1 > 0 :
            b[n - i - j - 2, :] = b[n - i - j -2, :] - b[n - j -1, :]*b[n - i - j - 2, n - j - 1]
            i += 1
            #print(b)
        j += 1
    #print(b)

    for i in range(n):
        c[i,0] = b[i, -1]
    return c
def L_decomp(d) :
    """
    :param d:   takes in  'N x N' or 'N x N+1' matrix 'd'
    :return:    lower triangular matrix(decomposed)
    """
    from numpy import identity, delete
    b = d.copy()
    c = identity(len(b[:, 0]))  # initializing an identity matrix
    n = len(b[:, 0])
    i = 0
    while i < n:
        j = 0
        while j + i + 1 < n:
            c[j + i + 1, i] = b[j + i + 1, i] / b[i, i]
            b[j + i + 1, :] = b[j + i + 1, :] - b[i, :] * float(b[j + i + 1, i] / b[i, i])
            j += 1
            # print(array(c))
            # print(b)
        i += 1

    #b = delete(b, n, 1)  # deleting the last column of the upper triangular matrix b
    return c
def U_decomp(d) :
    """
    :param d:   takes in 'N x N' or 'N x N+1' matrix 'd'
    :return:    'N x N' upper triangular matrix(decomposed)
    """
    from numpy import identity, delete
    b = d.copy()
    n = len(b[:, 0])
    i = 0
    while i < n:
        j = 0
        while j + i + 1 < n:
            b[j + i + 1, :] = b[j + i + 1, :] - b[i, :]*(b[j + i + 1, i] / b[i, i])
            j += 1
        i += 1
    if len(b[0,:]) + 1 == len(b[:,0]) :
        b = delete(b, n, 1)  # deleting the last column of the upper triangular matrix b
    return b
def matrix_inverse(b) :
    import sympy
    import numpy
    if determinant(b) == 0 :
        print("matrix inverse can't be calculated")
    else :
        d = b.copy()
        u = U_decomp(d)
        l = L_decomp(d)
        n = len(b[:,0])
        c = sympy.zeros(n,1)
        l = numpy.hstack((l,c))
        q = sympy.zeros(n,n)
        k = 0

        while k < n :
            l[k, -1] = 1
            w = gauss_el_roots(l)
            l[k,-1] = 0
            u = numpy.hstack((u,w))
            e = gauss_el_roots(u)
            i = 0
            while i < n :
                q[i,k] = e[i]
                i += 1
                #print(q)
            u = numpy.delete(u,n,1)
            k += 1
        return numpy.array(q)
def lagrange_interpolating_polynomial(a, x) :
    """
    :param a: a matrix with row elements as [x , y]
    :param x: point at which the function value needs to be calculated
    :return: the interpollated polynomial value at x
    """
    from sympy import ones
    try :
        i = 0
        li = ones(1, len(a[:, 0]))
        y = 0
        while i < len(a[:, 0]) :
            j = 0
            for j in range(len(a[:, 0])) :
                if i != j :
                    li[0, i] = (x - a[j, 0])/(a[i, 0] - a[j, 0])*li[0, i]
            y = y + li[0, i]*a[i, 1]
            i += 1
        return y
    except :
        "Something went wrong. Please note \n  1) 'a' is a matrix with row elements as [x , y] \n  2) x is the value at which the value needs to be calculated."
def lenier_spline(a,x) :

    try :
        i = 0
        while x > a[i,0] :
            i += 1
        y = ((a[i,1]-a[i-1,1])/(a[i,0]-a[i-1,0]))*(x - a[i,0]) +a[i,1]
    except :
        return "Something went wrong. Please note \n  1) 'a' is a matrix with row elements as [x , y] \n  2) x is the value at which the value needs to be calculated. \n  3) value of 'x' needs to be in the range of highest an lowest x value in the matrix"
    return y
def polynomial_regression(xx,yy,order) :
    """
    :param xx:      x values of the observation
    :param yy:      y values of the observation
    :param order:   order of the polynomial
    :return:        the polynmial equation of the order given
    """
    o = order
    c = array(zeros(o+1,o+2))
    i = 0
    while i < o+1 :
        j = 0
        while j < o+1 :
            c[i,j] = (sum((xx)**(i+j)))
            j += 1
        c[i,-1] = sum((xx**i)*yy)
        i += 1
    d = gauss_el_roots(c)
    x = symbols('x')
    i = 0
    y = 0
    while i < o + 1 :
        y = (d[i].item())*x**(i) + y
        i += 1
    return y
def determinant(a) :
    """
    :param a:   'N x N' matrix
    :return:    determinant of a
    """
    if len(a[:, 0]) != len(a[0, :]):
        print("\ngiven matrix is not an 'N x N' matrix")
    elif len(a[:,0]) == 2 == len(a[0,:]) :
        return a[0,0]*a[1,1] - a[1,0]*a[0,1]
    else:
        b = U_decomp(a)
        c = 1
        i = 0
        while i < len(b[0,:]) :
            c = c*b[i,i]
            i += 1
        return c
def hermite(y, n) :
    """
    :param y: point at which value of hermite polynomial needs to be calculated
    :param n: order of hermite polynomial
    :return: hermite polynomial
    """
    from sympy import symbols
    import sympy
    x = symbols('x')
    z = ((-1)**n)*(sympy.exp(x**2))*(sympy.diff(sympy.exp(-x**2), x, n))
    return z.subs(x,y)
def integration_trapezoidal(f,a,b,n) :
    """
    :param f: function of 'x'
    :param a: initial point
    :param b: final point
    :param n: number of parts
    :return: integral value
    """
    if n == 0 :
        print("n value cannot be zero")
    else :
        h = (b-a)/n
        m = 0
        i = h*(f(a) + f(b))/2
        while m < n :
            i = i + h*f(a + m*h)
            m += 1
        return i
def quadratic_spline(a,x0) :
    """
    :param a: matrix a representing values of x and y.
    :param x: point at which value is being calculated.
    :return: value of interpolated polynomial at x.
    """
    import sympy as sp
    import numpy as np
    i = 1
    y = ((a[1,1]-a[0,1])/(a[1,0]-a[0,0]))*(x0 - a[0,0]) + a[0,1]
    check = isinstance(x0,sp.symbol.Symbol)
    if check == True :
        return "work postponed until further notice"
    elif x0 < a[0,0] or x0 > a[-1,0] :
        return "x is out of bounds"
    elif x0 > a[0,0] and x0 < a[1,0] :
        return y
    else :
        x = sp.symbols('x')
        while x0 > a[i,0] :
            y = ((a[i, 1] - a[i-1, 1]) / (a[i, 0] - a[i-1, 0])) * (x - a[i-1, 0]) + a[i-1, 1]
            b = np.array([[2*a[i,0],1,0,sp.diff(y,x).subs(x,a[i,0])],[a[i,0]**2,a[i,0],1,a[i,1]],[a[i+1,0]**2,a[i+1,0],1,a[i+1,1]]])
            c = gauss_el_roots(b)
            i += 1
        return c[0,0]*x0**2 + c[1,0]*x0 + c[2,0]
def simpsons13_function(f,a,b,n) :
    """
        :param f: function of 'x'
        :param a: initial point
        :param b: final point
        :param n: number of parts
        :return: integral value
        """
    if n == 0 :
        print("n value cannot be zero")
    else :
        h = (b-a)/n
        m = 0
        i = h*(f(a) + f(b))/3
        while 2*m+1 < n :
            i = i + 4*h*f(a + (2*m+1)*h)/3
            m += 1
        m = 1
        while 2*m < n :
            i = i + 2*h*f(a + (2*m)*h)/3
            m += 1
        return i
def simpson38_function(f,a,b,n) :
    """
        :param f: function of 'x'
        :param a: initial point
        :param b: final point
        :param n: number of parts
        :return: integral value
        """
    if n <= 0 :
        print("n value cannot be zero or negative")
    elif n%3 == 0 :
        h = (b-a)/n
        m = 1
        i = 3*h*(f(a) + f(b))/8
        while 3*m < n :
            i = i + 3*h*f(a + (3*m)*h)/4
            m += 1
        m = 0
        while 3*m+1 < n :
            i = i + 9*h*f(a + (3*m+1)*h)/8
            m += 1
        m = 0
        while 3 * m + 2 < n:
            i = i + 9 * h * f(a + (3 * m + 2) * h) / 8
            m += 1
        return i
    elif n%3 == 2 :
        i = simpsons13_function(f,a,a+2*((b-a)/n),2)
        a = a+2*((b-a)/n)
        n = n - 2
        h = (b - a) / n
        m = 1
        i = i + 3 * h * (f(a) + f(b)) / 8
        while 3 * m < n:
            i = i + 3 * h * f(a + (3 * m) * h) / 4
            m += 1
        m = 0
        while 3 * m + 1 < n:
            i = i + 9 * h * f(a + (3 * m + 1) * h) / 8
            m += 1
        m = 0
        while 3 * m + 2 < n:
            i = i + 9 * h * f(a + (3 * m + 2) * h) / 8
            m += 1
        return i
    elif n%3 == 1:
        i = integration_trapezoidal(f,a,a+(b-a)/n,1)
        a = a+(b-a)/n
        n = n - 1
        h = (b - a) / n
        m = 1
        i = i + 3 * h * (f(a) + f(b)) / 8
        while 3 * m < n:
            i = i + 3 * h * f(a + (3 * m) * h) / 4
            m += 1
        m = 0
        while 3 * m + 1 < n:
            i = i + 9 * h * f(a + (3 * m + 1) * h) / 8
            m += 1
        m = 0
        while 3 * m + 2 < n:
            i = i + 9 * h * f(a + (3 * m + 2) * h) / 8
            m += 1
        return i

def euler(f,x0,tf,nMax) :
    """
    :param f: array of functions. argument of every element of f is an array of same sequence as x0.
    :param x0: initial values of arguments on which th euler method is to be applied.
    :param tf: final value of first argument of x0.
    :param nMax: number of data points to be made.
    :return: set of data points corresponding to order of initial values in x0.
    """
    import numpy as np
    h = (tf - x0[0])/nMax
    datalist = np.array([x0])
    while len(datalist[:,0]) < nMax :
        prev = datalist[-1,:]
        rate = np.array([])
        for func in f:
            rate = np.append(rate,func(prev))
        next = prev + h*(rate)
        datalist = np.append(datalist,[next],axis=0)
    return datalist


def euler_improved(f,x0,tf,nMax) :
    """
    :param f: array of functions. argument of every element of f is an array of same sequence as x0.
    :param x0: initial values of arguments on which th improved euler method is to be applied.
    :param tf: final value of first argument of x0.
    :param nMax: number of data points to be made.
    :return: set of data points corresponding to order of initial values in x0.
    """
    import numpy as np
    h = (tf - x0[0])/nMax
    datalist = np.array([x0])
    while len(datalist[:,0]) < nMax :
        prev = datalist[-1,:]
        rate1,rate2 = np.array([]),np.array([])
        for func in f:
            rate1 = np.append(rate1,func(prev))
        for func in f:
            rate2 = np.append(rate2,func(prev + h*rate1))
        next = prev + h*(rate1 + rate2)/2
        datalist = np.append(datalist,[next],axis=0)
    return datalist

def rk4_at_x(f,y0,x0,x,h) :
    """
    :param f: differential equation of derivative of y.(Takes in two arguements x and y)
    :param y0: value of x0 at y0
    :param x: value at which y needs to be found
    :param h:step height
    :return:value of y at x
    """
    n = int((x-x0)/h)
    y = y0
    i = 0
    while i < n:
        k1 = f(x0,y)
        k2 = f(x0 + h/2,y + k1*h/2)
        k3 = f(x0 + h/2,y + k2*h/2)
        k4 = f(x0 + h,y + k3*h)
        i += 1
        y = y + (k1 + 2*k2 + 2*k3 + k4)*h/6
        x0 = x0 + h
    return y
def rk4(f,x0,tf,nMax) :
    """
    :param f: array of functions. argument of every element of f is an array of same sequence as x0.
    :param x0: initial values of arguments on which th rk4 method is to be applied.
    :param tf: final value of first argument of x0.
    :param nMax: number of data points to be made.
    :return: set of data points corresponding to order of initial values in x0.
    """
    import numpy as np
    h = (tf - x0[0])/nMax
    datalist = np.array([x0])
    while len(datalist[:,0]) < nMax :
        prev = datalist[-1,:]
        rate1,rate2,rate3,rate4 = np.array([]),np.array([]),np.array([]),np.array([])
        for func in f:
            rate1 = np.append(rate1,func(prev))
        for func in f:
            rate2 = np.append(rate2,func(prev + h*rate1/2))
        for func in f:
            rate3 = np.append(rate3,func(prev + h*rate2/2))
        for func in f:
            rate4 = np.append(rate4,func(prev + h*rate3))
        next = prev + h*(rate1 + 2*(rate2 + rate3) + rate4)/6
        datalist = np.append(datalist,[next],axis=0)
    return datalist
