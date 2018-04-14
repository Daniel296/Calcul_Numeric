import numpy as np
import matplotlib.pyplot as plt

def readFromFile(file_path):
    fd = open(file_path, 'rb')
    values = {}
    """ Read the values of the function """
    while True:
        line = fd.readline()
        array = line.split()
        if len(array) < 2:
            break
        values[float(array[0])] = float(array[1])
    """ Get value of n """
    line = fd.readline().strip()
    n_array = line.split("=")
    if len(n_array) < 2:
        return False, None, None, None, None, None
    n = int(n_array[1])
    """ Get value of x0 """
    line = fd.readline().strip()
    x0_array = line.split("=")
    if len(x0_array) < 2:
        return False, None, None, None, None, None
    x0 = float(x0_array[1])
    """ Get value of xn """
    line = fd.readline().strip()
    xn_array = line.split("=")
    if len(xn_array) < 2:
        return False, None, None, None, None, None
    xn = float(xn_array[1])
    """ Get value of x """
    line = fd.readline().strip()
    x_array = line.split("=")
    if len(x_array) < 2:
        return False, None, None, None, None, None
    x = float(x_array[1])

    return True, values, n, x0, xn, x

def progresiveNewton(values, n, x0, xn, x):
    if xn < x0:
        return False, None
    h = (xn - x0)/(n-1)
    """ Determine the interpolation points """
    x_interpol = [0] * (n)
    x_interpol[0] = x0
    for i in range(1, n):
        x_interpol[i] = x0 + i * h
    """ Determine t """
    t = (x - x0)/h
    """ Start computing L """
    y = [0] * (n)
    y = [0] * (n)
    s = [0] * (n)
    y[0] = values[x0]
    s[0] = 1
    s[1] = t
    L = 0

    """ First step """
    for counter in range(1, n):
        y[counter] = values[values.keys()[counter + 1]] - values[values.keys()[counter]]

    """ Following steps """
    for pas in range(2, n):
        s[pas] = s[pas - 1] * (t - pas + 1) / pas

        z = [0] * (n-pas)
        for counter in range(0, n-pas):
            z[counter] = y[pas + counter] - y[pas + counter - 1]
        for counter in range(0, n-pas):
            y[pas + counter] = z[counter]

    """ Compute L """
    for counter in range(0, n):
        L += s[counter] * y[counter]

    return True, L

def interpolSmallestSquares(values, n, x0, xn, x):
    result = 0
    if xn < x0:
        return False, None
    h = (xn - x0) / (n - 1)
    """ Determine the interpolation points """
    x_interpol = [0] * (n)
    x_interpol[0] = x0
    for i in range(1, n):
        x_interpol[i] = x0 + i * h
    #Construct the polynom: I will have something like: S(x) = a1*pow(x, 4) + a2*pow(x,3) + a3*pow(x,2) + a4*pow(x,1) + a5
    #Construct the coefficient matrix
    #For each of the interpolation values
    coeff_matrix = np.empty([n, n])
    result_array = np.empty(n)
    for line_index in range(0, 5):
        for col_index in range(0, 5):
            coeff_matrix[line_index][col_index] =  pow(x_interpol[line_index], 4-col_index)
        result_array[line_index] = values[x_interpol[line_index]]

    #Determine the a coefficient array
    a = np.linalg.solve(coeff_matrix, result_array)
    #Determine the value in the point based on Horner's scheme
    d = np.empty(n)
    d[0] = a[0]
    for i in range(1, n):
        d[i] = a[i] + d[i-1] * x

    result = d[n-1]
    return True, result

def func(x, values):
    val = []
    for x in x:
        val.append(values[x])
    return val

def newton(x, values, n, x0, xn):
    val = []
    for x in x:
        val.append(progresiveNewton(values, n, x0, xn, x)[1])
    return val

def leastSquares(x, values, n, x0, xn):
    val = []
    for x in x:
        val.append(interpolSmallestSquares(values, n, x0, xn, x)[1])
    return val

def main():
    check, values, n, x0, xn, x = readFromFile("input.txt")
    if check == True:
        check, result = progresiveNewton(values, n, x0, xn, x)
        if check == True:
            print("Result with progressive Newton: {0}".format(result))
        check, result = interpolSmallestSquares(values, n, x0, xn, x)
        if check == True:
            print("Result with least squares interpolation: {0}".format(result))
    else:
        print("Couldn't find result!")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    precision = 50
    x = [1, 2, 3, 4, 5]
    ax.plot(x, func(x, values), color='g', label='Original', alpha=.3)
    x = np.arange(1, 100, 0.5)
    ax.plot(x, newton(x, values, n, x0, xn), color='b', label='Progressive Newton', alpha=.3)
    ax.plot(x, leastSquares(x, values, n, x0, xn), color='r', label='Least Squares Interpol', alpha=.3)
    ax.legend(fancybox=True)
    plt.show()

if __name__ == '__main__':
    main()
