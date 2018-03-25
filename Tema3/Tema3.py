EPS = pow(10, -7)

def getSize(vector):
    size_vector = len(vector)
    return size_vector;

def readFromFile(file_path, limit_line):
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
                    #print("I am on row " + str(line_index) + "col: " + str(column_index) + "val: " + str(value))
                    #print("Exists on column " + str(pair[1]) + "line: " + str(line_index) + "value: " + str(pair[0]))
                    already_exist = True
                    break
        if already_exist == True:
            #print("Added " + str(value) + " to " + str(pair[0]))
            pair[0] += value
        else:
            """ If I already have 10 not null elements on that line """
            if len(A[line_index]) == limit_line:
                print("Dimension of matrix exceeded!\n")
                exit(1)
            """ Else """
            A[line_index].append([value, column_index])

    return matrix_size, A, b

def isEqualVV(vector_1, vector_2):
    v_size = len(vector_1)
    if v_size != len(vector_2):
        return False
    else:
        for list_index in range(v_size):
            if vector_1[list_index] - vector_2[list_index] > EPS:
                return False
    return True

def isEqualMM(matrix_1, matrix_2):
    matrix_size = len(matrix_1)

    for line_index in range(matrix_size):
        if getSize(matrix_1[line_index]) != getSize(matrix_2[line_index]):
            return False
        for pair in matrix_1[line_index]:
            value_1 = pair[0]
            col_1 = pair[1]
            found_col_1 = False
            for pair_2 in matrix_2[line_index]:
                if col_1 == pair_2[1]:
                    found_col_1 = True
                    if pair_2[0] - value_1 > EPS:
                        return False
            if found_col_1 == False:
                return False
    return True

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

def sumMatrixMatrix(A, B):
    matrix_size = len(A)
    """ Initialize the result matrix """
    result = [[] for _ in range(matrix_size)]
    """ For each line in the matrix """
    for line_index in range(matrix_size):
        size_A = getSize(A[line_index])
        size_B = getSize(B[line_index])
        d = {}
        for el_index in range(size_A):
            """ If already have value, add the new value to the old one """
            if d.get(A[line_index][el_index][1]):
                d[A[line_index][el_index][1]] += A[line_index][el_index][0]
            else:
                d[A[line_index][el_index][1]] = A[line_index][el_index][0]
        for el_index in range(size_B):
            if d.get(B[line_index][el_index][1]):
                d[B[line_index][el_index][1]] += B[line_index][el_index][0]
            else:
                d[B[line_index][el_index][1]] = B[line_index][el_index][0]

        for j in sorted(d.keys()):
            result[line_index].append((d[j], j))

    return result

def main():
    matrix_size, A, b = readFromFile('a.txt', 10)
    x = [matrix_size - i for i in range(2018)]
    result = mulMatrixVector(A, x)
    """ Check that Matrix*Vector multiplication is equal to b vector """
    print(isEqualVV(result, b))
    """ Check that A_a + A_b is equal to matrix from aplusb.txt """
    matrix_size_b, A_b, b_b = readFromFile('b.txt', 10)
    AplusB_size, AplusB_matrix, AplusB_vector = readFromFile('aplusb.txt', 20)
    result = sumMatrixMatrix(A, A_b)
    print(isEqualMM(result, AplusB_matrix))

if __name__ == '__main__':
    main()
