import numpy as np
from copy import deepcopy
from scipy.spatial import distance
from math import sqrt

""" Aux functions """


def vectorXmatrix(vector, matrix):
    if (len(vector) != len(matrix)):
        print("ERROR nr elements vector != nr lines matrix!")
        return []

    result = []
    for nrCol in range(0, len(matrix[0])):
        element = 0.0
        for nrLine in range(0, len(matrix)):
            element += vector[nrLine] * matrix[nrLine][nrCol]
        result.append(element)

    return result


def matrixXvector(matrix, vector):
    if (len(matrix[0]) != len(vector)):
        print("ERROR nr elements vector != nr columns matrix!")
        return []

    result = []

    for nrLine in range(0, len(matrix)):
        element = 0.0
        for nrCol in range(0, len(matrix[0])):
            element += vector[nrCol] * matrix[nrLine][nrCol]
        result.append(element)

    return result


def euclidianNorma(z):
    norma = 0.0
    for i in range(0, len(z)):
        norma = norma + z[i] * z[i]
    return sqrt(norma)


def vectorMinusVector(a, b):
    result = []
    for i in range(0, len(a)):
        result.append(a[i] - b[i])
    return result


""" Bonus - Thomas Algorithm """


def swap(a, b, poz_a, poz_b):
    aux = a[poz_a]
    a[poz_a] = b[poz_b]
    b[poz_b] = aux
    return a, b


def tri_diag_matrix_solver(a, b, c, results):
    diag_size = len(a)  # number of equations
    d, e, cc, results_copy = map(np.array, (a, b, c, results))  # copy arrays
    f = np.zeros(diag_size - 1)
    for it in range(0, diag_size - 1):
        if d[it] < cc[it]:
            d, cc = swap(d, cc, it, it)
            if it < diag_size - 1:
                d, e = swap(d, e, it + 1, it)
            if it < diag_size - 2:
                e, f = swap(e, f, it + 1, it)
            aux = float(results_copy[it])
            results_copy[it] = float(results_copy[it + 1])
            results_copy[it + 1] = aux

        temp = (-1 * d[it] / cc[it])
        d[it + 1] = temp * d[it + 1] + e[it]
        if it < diag_size - 2:
            e[it + 1] = temp * e[it + 1] + f[it]
        if it < diag_size - 3:
            f[it + 1] = f[it + 1] * temp
        results_copy[it + 1] = results_copy[it] + results_copy[it + 1] * temp

    xc = np.zeros(diag_size)
    xc[diag_size - 1] = results_copy[diag_size - 1] / d[diag_size - 1]
    xc[diag_size - 2] = (results_copy[diag_size - 2] - xc[diag_size - 1] * e[diag_size - 2]) / d[diag_size - 2]
    for il in range(diag_size - 3, -1, -1):
        xc[il] = (results_copy[il] - xc[il + 1] * e[il] - xc[il + 2] * f[il]) / d[il]

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


# def pivot_test_bonus(a, b, c):


if __name__ == "__main__":
    fd = open('input', 'r')
    buf = fd.readlines()
    fd.close()
    n = int(buf[0].strip())
    epsilon = pow(10, -int(buf[1].strip()))

    A = np.zeros((n, n))
    for line in range(2, n + 2):
        A[line - 2] = np.array([float(i.strip()) for i in buf[line].split(' ')])
    b = np.array([float(i.strip()) for i in buf[-1].split(' ')])
    Acopy = deepcopy(A)
    Acopy_2 = deepcopy(A)
    bcopy = deepcopy(b)
    bcopy_2 = deepcopy(b)

    my_result = partial_gauss(n, epsilon, A, b)
    np_result = (np.linalg.solve(Acopy, bcopy))
    print("Rezultatul meu:{0}".format(my_result))
    print("Rezultatul cu numpy:{0}".format(np_result))
    result_mul = matrixXvector(Acopy, my_result)
    result_sub = vectorMinusVector(result_mul, bcopy)
    norma = euclidianNorma(result_sub)
    # euclidean_distance = distance.euclidean(np.dot(Acopy, my_result), bcopy)
    print("Verificarea:{0}".format(norma))
    inverse = np.linalg.inv(Acopy)
    print("Inversa cu numpy:{0}".format(inverse))
    result_sub = vectorMinusVector(my_result, np_result)
    norma = euclidianNorma(result_sub)
    # euclidean_distance = distance.euclidean(my_result, np_result)
    print("Norma euclidiana intre rezultate:{0}".format(norma))
    result_mul = matrixXvector(inverse, bcopy)
    result_sub = vectorMinusVector(my_result, result_mul)
    norma = euclidianNorma(result_sub)
    # euclidean_distance = distance.euclidean(my_result, np.dot(inverse, bcopy))
    print("Norma euclidiana intre my_result si rezultatul asteptat:{0}".format(norma))

    # For tridiag matrix
    a = [Acopy_2[i][i] for i in range(n)]
    b = [Acopy_2[i][i + 1] for i in range(n - 1)]
    c = [Acopy_2[i + 1][i] for i in range(n - 1)]
    d = deepcopy(bcopy_2)
    print 'Tridiag matrix solution: ', tri_diag_matrix_solver(a, b, c, d)
