import numpy as np
from copy import deepcopy
from scipy.spatial import distance

def cauta_pivot(l, epsilon, A, b):
    maxindex = abs(A[l:, l]).argmax() + l
    if abs(A[maxindex, l]) <= epsilon:
        raise ValueError("Matrix is singular.")
    # Swap rows
    if maxindex != l:
        A[[l, maxindex]] = A[[maxindex, l]]
        b[[l, maxindex]] = b[[maxindex, l]]
    return A, b

def Gauss_Partial(n, epsilon, A, b):
    '''
    Gaussian elimination with partial pivoting.
    % input: n is the dimension of the matrix
    %        A is an n x n nonsingular matrix
    %        b is an n x 1 vector
    %        epsilon is the precision of the computations
    % output: x is the solution of Ax=b.
    '''
    l = 0
    for l in range(n-1):
        #Choose pivot
        A, b = cauta_pivot(l, epsilon, A, b)
        for row in range(l+1, n):
            multiplier = A[row][l]/A[l][l]
            #the only one in this column since the rest are zero
            A[row][l] = multiplier
            for col in xrange(l + 1, n):
                A[row][col] = A[row][col] - multiplier*A[l][col]
            #Equation solution column
            b[row] = b[row] - multiplier*b[l]
    x = np.zeros(n)
    l = n-1
    x[l] = b[l]/A[l,l]
    while l >= 0:
        x[l] = (b[l] - np.dot(A[l,l+1:],x[l+1:]))/A[l,l]
        l = l-1
    return x

if __name__ == "__main__":
    A = np.array([[1., -1., 1., -1.], [1., 0., 0., 0.], [1., 1., 1., 1.], [1., 2., 4., 8.]])
    Acopy = deepcopy(A)
    b = np.array([[14.], [4.], [2.], [2.]])
    bcopy = deepcopy(b)
    epsilon = pow(10, -10);
    n = len(A)
    my_result = Gauss_Partial(n, epsilon, A, b)
    np_result = (np.linalg.solve(Acopy, bcopy))
    print("Rezultatul meu:{0}".format(my_result))
    print("Rezultatul cu numpy:{0}".format(np_result))
    euclidean_distance = distance.euclidean(np.dot(Acopy, my_result), bcopy)
    print("Verificarea:{0}".format(euclidean_distance))
    inverse = np.linalg.inv(Acopy)
    print("Inversa cu numpy:{0}".format(inverse))
    euclidean_distance = distance.euclidean(my_result, np_result)
    print("Norma euclidiana intre rezultate:{0}".format(euclidean_distance))
    euclidean_distance = distance.euclidean(my_result, np.dot(inverse, bcopy))
    print("Norma euclidiana intre my_result si rezultatul asteptat:{0}".format(euclidean_distance))


