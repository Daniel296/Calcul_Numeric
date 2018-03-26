from math import sqrt

EPS = pow(10, -7)
MAX_K = 10000
MAX_DELTA = pow(10, 8)

def euclideanNorm(Z):
    norm = 0.0
    for index_vector in range(0, len(Z)):
        norm = norm + Z[index_vector] * Z[index_vector]
    return sqrt(norm)

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

def getValue(A, line, column):
    for index_vector in range(len((A[line]))):
        if A[line][index_vector][1] == column:
            return A[line][index_vector][0];
    return -200000

def getNewX(A, b, xGS):
    delta = 0
    for line_index in range(len(b)):
        sum_bellow_diagonal = 0
        sum_above_diagonal = 0
        for col_index in range(line_index):
            temp_val = getValue(A, line_index, col_index)
            if temp_val != -200000:
                sum_bellow_diagonal += temp_val * xGS[col_index]
        for col_index in range(line_index + 1, len(b)):
            temp_val = getValue(A, line_index, col_index)
            if temp_val != -200000:
                sum_above_diagonal += temp_val * xGS[col_index]
        delta += pow(((b[line_index] - sum_bellow_diagonal - sum_above_diagonal)/getValue(A, line_index, line_index)) - xGS[line_index], 2)
        xGS[line_index] = (b[line_index] - sum_bellow_diagonal - sum_above_diagonal)/getValue(A, line_index, line_index)
    return sqrt(delta), xGS

def Gauss_Siedel(A, b):
    xGS = []
    k = 0
    for index_vector in range(len(b)):
        xGS.append(0.0)
    while True:
        delta, xGS = getNewX(A, b, xGS)
        k += 1
        print("Delta: {0}".format(delta))
        print("k: {0}".format(k))
        if (delta < EPS or k > MAX_K or delta > MAX_DELTA):
            break
    if delta < EPS:
        return True, xGS
    else:
        return False

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
                if column_index == line_index:
                    if abs(value) < EPS:
                        print("Value 0 found on diagonal, for column index equal to: {0} and line index equal to: {1}".format(line_index, column_index))
                        exit(0)
        if already_exist == True:
            pair[0] += value
        else:
            A[line_index].append([value, column_index])

    return matrix_size, A, b

def main():
    matrix_size, A, b = readFromFile("m_rar_2018_5.txt")
    convergent, xGS = Gauss_Siedel(A, b)
    if (convergent):
        print("Convergent!\n")
        result = mulMatrixVector(A, xGS)
        print(euclideanNorm(subMatrices(result, b)))
    else:
        print("Divergent!\n")

if __name__ == '__main__':
    main()
