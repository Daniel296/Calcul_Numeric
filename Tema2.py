import numpy as np
from copy import deepcopy
from scipy.spatial import distance
from math import sqrt

""" Aux functions """
def vectorXmatrix(vector, matrix):
    if (len(vector) != len(matrix)):
        print("ERROR nr elements vector != nr lines matrix!");
        return [];

    result = [];
    for nrCol in range(0, len(matrix[0])):
        element = 0.0;
        for nrLine in range(0, len(matrix)):
            element += vector[nrLine] * matrix[nrLine][nrCol];
        result.append(element);

    return result;

def matrixXvector(matrix, vector):
    if (len(matrix[0]) != len(vector)):
        print("ERROR nr elements vector != nr columns matrix!");
        return [];

    result = [];

    for nrLine in range(0, len(matrix)):
        element = 0.0;
        for nrCol in range(0, len(matrix[0])):
            element += vector[nrCol] * matrix[nrLine][nrCol];
        result.append(element);

    return result;

def euclidianNorma(z):
    norma = 0.0;
    for i in range(0, len(z)):
        norma = norma + z[i] * z[i];
    return sqrt(norma);

def vectorMinusVector(a, b):
    result = [];
    for i in range(0, len(a)):
        result.append(a[i] - b[i]);
    return result;

""" Bonus - Thomas Algorithm """
def tri_diag_matrix_solver(a, b, c, results):
    nf = len(results)  # number of equations
    ac, bc, cc, results_copy = map(np.array, (a, b, c, results))  # copy arrays
    for it in xrange(1, nf):
        mc = ac[it - 1] / bc[it - 1]
        bc[it] = bc[it] - mc * cc[it - 1]
        results_copy[it] = results_copy[it] - mc * results_copy[it - 1]

    xc = bc
    xc[-1] = results_copy[-1] / bc[-1]

    for il in xrange(nf - 2, -1, -1):
        xc[il] = (results_copy[il] - cc[il] * xc[il + 1]) / bc[il]

    return xc

""" Function for pivot searching """
def search_pivot(l, epsilon, A, b):
    maxindex = abs(A[l:, l]).argmax() + l
    if abs(A[maxindex, l]) <= epsilon:
        raise ValueError("Matrix is singular.")
    # Swap rows
    if maxindex != l:
        A[[l, maxindex]] = A[[maxindex, l]]
        b[[l, maxindex]] = b[[maxindex, l]]
    return A, b


""" Gauss cu pivotare partiala """
def partial_gauss(n, epsilon, A, b):
    for l in range(n - 1):
        # Choose pivot
        A, b = search_pivot(l, epsilon, A, b)
        for row in range(l + 1, n):
            if abs(A[l][l]) > epsilon:
                ratia = A[row][l] / A[l][l]
                # the only one in this column since the rest are zero
                A[row][l] = ratia
                for col in range(l + 1, n):
                    A[row][col] = A[row][col] - ratia * A[l][col]
                # Equation solution column
                b[row] = b[row] - ratia * b[l]
            else:
                print('Division by 0! Returning...')
                exit(1)
    x = np.zeros(n)
    l = n - 1
    x[l] = b[l] / A[l, l]
    while l >= 0:
        x[l] = (b[l] - np.dot(A[l, l + 1:], x[l + 1:])) / A[l, l]
        l = l - 1

    return x


if __name__ == "__main__":
    fd = open('input', 'r')
    buf = fd.readlines()
    fd.close()
    n = int(buf[0].strip())
    epsilon = pow(10, -int(buf[1].strip()))

    A= np.zeros((n, n))
    for line in range(2, n + 2):
        A[line - 2] = np.array([float(i.strip()) for i in buf[line].split(' ')])
    b = np.array([float(i.strip()) for i in buf[-1].split(' ')])
    Acopy = deepcopy(A)
    bcopy = deepcopy(b)

    my_result = partial_gauss(n, epsilon, A, b)
    np_result = (np.linalg.solve(Acopy, bcopy))
    print("Rezultatul meu:{0}".format(my_result))
    print("Rezultatul cu numpy:{0}".format(np_result))
    result_mul = matrixXvector(Acopy, my_result)
    result_sub = vectorMinusVector(result_mul, bcopy)
    norma = euclidianNorma(result_sub)
    #euclidean_distance = distance.euclidean(np.dot(Acopy, my_result), bcopy)
    print("Verificarea:{0}".format(norma))
    inverse = np.linalg.inv(Acopy)
    print("Inversa cu numpy:{0}".format(inverse))
    result_sub = vectorMinusVector(my_result, np_result)
    norma = euclidianNorma(result_sub)
    #euclidean_distance = distance.euclidean(my_result, np_result)
    print("Norma euclidiana intre rezultate:{0}".format(norma))
    result_mul = matrixXvector(inverse, bcopy)
    result_sub = vectorMinusVector(my_result, result_mul)
    norma = euclidianNorma(result_sub)
    #euclidean_distance = distance.euclidean(my_result, np.dot(inverse, bcopy))
    print("Norma euclidiana intre my_result si rezultatul asteptat:{0}".format(norma))

    a = [1, 2, 3]
    b = [1, 1, 2, 2]
    c = [7, 8, 9]
    d = [5, 5, 5, 5]
    #print(tri_diag_matrix_solver(a, b, c, d))
