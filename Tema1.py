import numpy as np

""" Problema 1 """
def machinePrecisionSum():
    m = -1
    while (True):
        if (1.0 + pow(10, m) == 1.0):
            return pow(10, m + 1), (m + 1)
        m -= 1;
#print("Final result is stated for u = {0} and power = {1}".format(machinePrecisionSum()[0], machinePrecisionSum()[1]))

""" Problema 2 """
def checkNonAssociativity(x = 1.0):
    y = z = machinePrecisionSum()[0]
    boolean_value = ((x + y) + z == x + (y + z))
    if boolean_value is False:
        print("The operation is non-associative")
    else:
        print("The operation is associative")
#print(checkNonAssociativity())

def checkMultiplicationAssociativity():
    y = z = machinePrecisionSum()[0]
    x = 1.0
    boolean_value = ((x * y) * z == x * (y * z))
    while boolean_value:
        x = x * 1.1
        boolean_value = ((x * y) * z == x * (y * z))
    return x, (x * y) * z, x * (y * z)
#print(checkMultiplicationAssociativity())

""" Problema 3 """
"""
    Implementation of the strassen algorithm, similar to
    http://en.wikipedia.org/w/index.php?title=Strassen_algorithm&oldid=498910018#Source_code_of_the_Strassen_algorithm_in_C_language
"""
def strassen(A, B):
    n = len(A)
    if n==1:
        C = np.multiply(A, B)
    else:
        new_size = n//2
        """ Initialize the submatrices """
        a11 = np.zeros(shape=(new_size, new_size))
        a12 = np.zeros(shape=(new_size, new_size))
        a21 = np.zeros(shape=(new_size, new_size))
        a22 = np.zeros(shape=(new_size, new_size))
        b11 = np.zeros(shape=(new_size, new_size))
        b12 = np.zeros(shape=(new_size, new_size))
        b21 = np.zeros(shape=(new_size, new_size))
        b22 = np.zeros(shape=(new_size, new_size))

        aResult = np.zeros(shape=(new_size, new_size))
        bResult = np.zeros(shape=(new_size, new_size))

        """ Devide the matrices in 4 submatrices """
        for i in range(0, new_size):
            for j in range(0, new_size):
                """ Bottom left """
                a11[i][j] = A[i][j]
                """ Bottom right """
                a12[i][j] = A[i][j + new_size]
                """ Down left """
                a21[i][j] = A[i + new_size][j]
                """ Down size """
                a22[i][j] = A[i + new_size][j + new_size]

                """ Bottom left """
                b11[i][j] = B[i][j]
                """ Bottom right """
                b12[i][j] = B[i][j + new_size]
                """ Down left """
                b21[i][j] = B[i + new_size][j]
                """ Down size """
                b22[i][j] = B[i + new_size][j + new_size]

                aResult = np.add(a11, a22)
                bResult = np.add(b11, b22)
                p1 = strassen(aResult, bResult)
                
                aResult = np.subtract(a21, a22)
                p2 = strassen(aResult, b11)
                
                bResult = np.subtract(b12, b22)
                p3 = strassen(a11, bResult)
                
                bResult = np.subtract(b21, b11)
                p4 = strassen(a22, bResult)
                
                aResult = np.add(a11, a12)
                p5 = strassen(aResult, b22)
                
                aResult = np.subtract(a21, a11)
                bResult = np.add(b11, b12)
                p6 = strassen(aResult, bResult)
                
                aResult = np.subtract(a12, a22)
                bResult = np.subtract(b21, a22)
                p7 = strassen(aResult, bResult)
                
                c11 = p1 + p4 - p5 + p7
                c12 = p3 + p5
                c21 = p2 + p4
                c22 = p1 + p3 - p2 + p6
                
                C = np.zeros(shape=(n, n))
                for i in range(0, new_size):
                    for j in range(0, new_size):
                        C[i][j] = c11[i][j]
                        C[i][j + new_size] = c12[i][j]
                        C[i + new_size][j] = c21[i][j]
                        C[i + new_size][j + new_size] = c22[i][j]
                        
                return C

if __name__ == "__main__":                
    A = np.zeros(shape=(8, 8))
    B = np.zeros(shape=(8, 8))
    for i in range(0, 8):
        for j in range(0, 8):
            A[i][j] = i
            B[i][j] = j          
    print(strassen(A, B))