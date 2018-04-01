import numpy as np

EPS = pow(10, -5)
MAX_K = 10000
MAX_DELTA = pow(10, 8)

def subMatrices(A, B):
    result = [];
    for index_vector in range(len(A)):
        result.append(A[index_vector] - B[index_vector]);
    return result;

def mulMatrixVector(A, vector):
    matrix_size = len(A)
    """ Initialize the result matrix """
    result = [0 for _ in range(matrix_size)]

    """ For each row, multiply the not null elements with the corresponding elements in the vector """
    for line_index in range(matrix_size):
        total_mul_line = 0
        for pair in A[line_index]:
            value = pair[0]
            column_index = pair[1]
            total_mul_line += value * vector[column_index]
        result[line_index] = total_mul_line

    return result

def verifyMatrix(A):
    for line in range(len(A)):
        foundInLine = False
        for column in range(len(A[line])):
            """ If the line is equal to the column """
            if A[line][column][1] == line:
                foundInLine = True
                if abs(A[line][column][0]) < EPS:
                    print("Can not apply Gauss_Siedel!\n")
                    return False
        if foundInLine == False:
            print("Can not apply Gauss_Siedel!\n")
            return False
    return True

def getValue(A, line, column):
    for index_vector in range(len((A[line]))):
        if A[line][index_vector][1] == column:
            return A[line][index_vector][0];
    return -200000

def Gauss_Siedel(matrice, b, matrix_size):
    xGS = []
    for index_line in range(matrix_size):
        xGS.append(0.0)

    """ For each iteration """
    for k in range(MAX_K):
        is_conv = 1
        is_div = 1
        """ For each line """
        for index_line in range(matrix_size):
            suma_bellow_diagonal = 0
            suma_above_diagonal = 0
            contor = 0
            """ Get the two sums """
            while contor < len(matrice[index_line]) and  matrice[index_line][contor][1] < index_line:
                suma_bellow_diagonal += matrice[index_line][contor][0] * xGS[matrice[index_line][contor][1]]
                contor += 1
            contor += 1
            while contor < len(matrice[index_line]) and matrice[index_line][contor][1] < matrix_size:
                suma_above_diagonal += matrice[index_line][contor][0] * xGS[matrice[index_line][contor][1]]
                contor += 1
            """ Get the value of the element """
            elem = (b[index_line] - suma_bellow_diagonal - suma_above_diagonal) / getValue(matrice, index_line, index_line)
            """ Check if the ending conditions are met """
            if abs(xGS[index_line] - elem) >= EPS:
                is_conv = 0
            if abs(xGS[index_line] - elem) <= pow(10, 8):
                print("Delta exceeded!\n")
                is_div = 0
            """ Update the value of the current element in the vector """
            xGS[index_line] = elem
        if is_conv or is_div:
            break
    if is_conv and k != MAX_K - 1:
        print("Iteratia:{0} ".format(k))
        print("Convergent!\n")
        return xGS
    else:
        print("Iteratia:{0} ".format(k))
        print("Divergent!\n")
        return None

def readFromFile(file_path):
    fd = open(file_path, 'rb')
    """ Get matrix size """
    matrix_size = int(fd.readline().strip())
    """ Skip the white line between matrix size and elements of vector b """
    fd.readline()
    """ Get elements of b """
    b = []
    for i in range(matrix_size):
        b.append(float(fd.readline().strip()))
    """ Skip the white line between elements of b and elements of a """
    fd.readline()
    """ Get elements of A """
    """ Initializarea matricii A"""
    A = [[] for _ in range(matrix_size)]
    while True:
        line = fd.readline()
        if line == '':
            break
        """ Get all the elements in one line """
        el_in_line = line.strip().split(', ')
        value = float(el_in_line[0])
        line_index = int(el_in_line[1])
        column_index = int(el_in_line[2])
        already_exist = False
        if len(A[line_index]) > 0:
            for pair in A[line_index]:
                """ If I have an element on the same position as another one, I add the two values """
                if column_index == pair[1]:
                    already_exist = True
                    break
        if already_exist == True:
            pair[0] += value
        else:
            if len(A[line_index]) == 10:
                print("Dimension of matrix exceeded!\n")
                exit(1)
            A[line_index].append([value, column_index])

    return matrix_size, A, b

def multiply(A, x):
    return_vector = []
    for i in range(len(x)):
        temp_result = 0
        for index_vector in range(len(A[i])):
            temp_result += x[A[i][index_vector][1]] * A[i][index_vector][0]
        return_vector.append(temp_result)

    return np.array(return_vector)

def conjugate_gradient_method(A, b):
    x = []
    b = np.array(b)
    for i in range(len(b)):
        x.append(0.0)
    x = [x]
    r = []
    p = []
    r.append(b - multiply(A, x[0]))
    p.append(b - multiply(A, x[0]))
    k = 0
    while k < MAX_K:
        print(k)
        Ap = multiply(A, p[k])
        alpha = r[k].transpose().dot(r[k])
        alpha /= p[k].transpose().dot(Ap)
        x.append(x[k] + alpha * p[k])
        r.append(r[k] - alpha * Ap)
        if abs(r[k + 1].transpose().dot(r[k + 1])) < EPS:
            return x[k + 1]
        beta = r[k + 1].transpose().dot(r[k + 1])
        beta /= r[k].transpose().dot(r[k])
        p.append(r[k + 1] + beta * p[k])
        k += 1
    return x[k]

def main():
    matrix_size, A, b = readFromFile(r"m_rar_2018_5.txt")
    x = conjugate_gradient_method(A, b)
    print(x)
    print(np.linalg.norm(multiply(A, x) - b, np.inf))

if __name__ == '__main__':
    main()
