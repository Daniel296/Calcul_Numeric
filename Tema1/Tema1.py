import numpy as np

""" Problema 1 """


def machine_precision_sum():
    m = -1
    while True:
        if 1.0 + pow(10, m) == 1.0:
            return pow(10, m + 1), (m + 1)
        m -= 1


""" Problema 2 """


def check_non_associativity():
    x = 1.0
    y = z = machine_precision_sum()[0]
    boolean_value = ((x + y) + z == x + (y + z))
    if boolean_value is False:
        return "The sum operation is non-associative"
    else:
        return "The sum operation is associative"


def check_multiplication_associativity():
    y = z = machine_precision_sum()[0]
    x = 1.0
    boolean_value = ((x * y) * z == x * (y * z))
    while boolean_value:
        x = x * 1.1
        boolean_value = ((x * y) * z == x * (y * z))
    return x, (x * y) * z, x * (y * z)


""" Problema 3 """
"""
    Implementation of the strassen algorithm, similar to
    http://en.wikipedia.org/w/index.php?title=Strassen_algorithm&oldid=498910018#Source_code_of_the_Strassen_algorithm_in_C_language
"""


def multiply_strassen(matrix_a, matrix_b, n, n_min):
    if n <= n_min:
        return np.dot(matrix_a, matrix_b)

    new_size = n // 2
    """ Initialize the submatrices """
    a11 = np.zeros(shape=(new_size, new_size))
    a12 = np.zeros(shape=(new_size, new_size))
    a21 = np.zeros(shape=(new_size, new_size))
    a22 = np.zeros(shape=(new_size, new_size))
    b11 = np.zeros(shape=(new_size, new_size))
    b12 = np.zeros(shape=(new_size, new_size))
    b21 = np.zeros(shape=(new_size, new_size))
    b22 = np.zeros(shape=(new_size, new_size))

    """ Devide the matrices in 4 submatrices """
    for i in range(0, new_size):
        for j in range(0, new_size):
            """ Bottom left """
            a11[i][j] = matrix_a[i][j]
            """ Bottom right """
            a12[i][j] = matrix_a[i][j + new_size]
            """ Down left """
            a21[i][j] = matrix_a[i + new_size][j]
            """ Down size """
            a22[i][j] = matrix_a[i + new_size][j + new_size]

            """ Bottom left """
            b11[i][j] = matrix_b[i][j]
            """ Bottom right """
            b12[i][j] = matrix_b[i][j + new_size]
            """ Down left """
            b21[i][j] = matrix_b[i + new_size][j]
            """ Down size """
            b22[i][j] = matrix_b[i + new_size][j + new_size]

    p1 = multiply_strassen((a11 + a22), (b11 + b22), new_size, n_min)
    p2 = multiply_strassen((a21 + a22), b11, new_size, n_min)
    p3 = multiply_strassen(a11, (b12 - b22), new_size, n_min)

    p4 = multiply_strassen(a22, (b21 - b11), new_size, n_min)
    p5 = multiply_strassen((a11 + a12), b22, new_size, n_min)
    p6 = multiply_strassen((a21 - a11), (b11 + b12), new_size, n_min)
    p7 = multiply_strassen((a12 - a22), (b21 + b22), new_size, n_min)

    c11 = p1 + p4 - p5 + p7
    c12 = p3 + p5
    c21 = p2 + p4
    c22 = p1 + p3 - p2 + p6

    matrix_c = np.zeros(shape=(n, n))
    for i in range(0, new_size):
        for j in range(0, new_size):
            matrix_c[i][j] = c11[i][j]
            matrix_c[i][j + new_size] = c12[i][j]
            matrix_c[i + new_size][j] = c21[i][j]
            matrix_c[i + new_size][j + new_size] = c22[i][j]

    return matrix_c


if __name__ == "__main__":

    print("Final result is stated for u = {0} and power = {1}".format(machine_precision_sum()[0],
                                                                      machine_precision_sum()[1]))
    print(check_non_associativity())
    print("The product operation is associative for : ", check_multiplication_associativity())

    file = open("test.txt", "r")
    text = file.read()
    file.close()
    matrices = text.split("\n\n")
    if len(matrices) > 2 or  len(matrices) < 2:
        print("Please provide the proper number of matrices!\n");
    else:
        n_min = 1
        for matrix_index in range(0, len(matrices)):
            lines = matrices[matrix_index].split('\n')
            lines[:] = [item for item in lines if item != '']
            if matrix_index == 0:
                A = np.zeros(shape=(len(lines), len(lines)))
            elif matrix_index == 1:
                B = np.zeros(shape=(len(lines), len(lines)))
            for line_index in range(0, len(lines)):
                cols = lines[line_index].split()
                for col_index in range(0, len(cols)):
                    if matrix_index == 0:
                        A[line_index][col_index] = cols[col_index]
                    elif matrix_index == 1:
                        B[line_index][col_index] = cols[col_index]

        print(multiply_strassen(A, B, len(lines), n_min))
