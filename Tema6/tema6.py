import numpy as np
from random import sample

epsilon = pow(10, -9)
kmax = 1000000

def generateRareMatrix(n, sparsity_level):
    sample_pool = list(range(n))
    delta = 2000
    A = [[] for _ in range(n)]
    """ For each line """
    for i in range(n-1, 0, -1):
        """ Randomly generate the columns on which we have non zero elements """
        """ At each step I provide fewer values for the sample_pool """
        cols = sample(sample_pool, int((1 - sparsity_level) * (i + 1)))
        for col in cols:
            val = np.random.rand() * 2000
            A[i].append([val, col])
            """ Construct the superior half of the matrix. If we are on the diagonal, we add just one time the value """
            if i!=col:
                A[col].append([val, i])
        """ Pop the last elements because on the previous line I will go with 1 step less """
        sample_pool.pop()
    return A

def writeToFileGeneratedMatrix(A, file_path):
    fd = open(file_path, 'w')
    fd.write(str(len(A)) + "\n" + "\n")
    for line_index in range(0, len(A)):
        line = A[line_index]
        for pair in line:
            val = pair[0]
            col = pair[1]
            temp_str = str(val) + ', ' + str(line_index) + ', ' + str(col) + "\n"
            fd.write(temp_str)
    fd.close()

def getValue(A, line, col):
    for pair in A[line]:
        if pair[1] == col:
            return pair[0]
    return None

def checkMatrixSymmetry(A):
    for line_index in range(0, len(A)):
        line = A[line_index]
        for pair in line:
            value = pair[0]
            col = pair[1]
            transpose_val = getValue(A, col, line_index)
            if transpose_val != None:
                if value - transpose_val > epsilon:
                    return False
    return True

def readFromFile(file_path):
    fd = open(file_path, 'rb')
    """ Get matrix size """
    matrix_size = int(fd.readline().strip())
    """ Skip the white line between matrix size and elements of A """
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
            A[line_index].append([value, column_index])

    return matrix_size, A

def readSVD(file_path):
    fd = open(file_path, 'rb')
    """ Get number of lines """
    lines = int(fd.readline().strip())
    """ Skip the white line """
    fd.readline()
    """ Get number of columns """
    cols = int(fd.readline().strip())
    """ Skip the white line """
    fd.readline()
    """ Initialize matrix A """
    A = np.zeros(shape=(lines, cols))
    """ Get elements of A """
    line_index = 0;
    while True:
        line = fd.readline()
        col_index = 0;
        if line == '':
            break
        """ Get elements on each row """
        col_els = line.strip().split(' ')
        for value in col_els:
            A[line_index][col_index] = value
            col_index += 1
        line_index += 1

    return A

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

def get_eigenvalues(matrix_size, A):
    v = np.random.rand(matrix_size)
    v = np.dot((1/np.linalg.norm(v, ord=2)), v)
    w = mulMatrixVector(A, v)
    lda = np.inner(w, v)
    k = 0
    while k < kmax:
        v = np.dot((1/np.linalg.norm(w, ord=2)), w)
        w = mulMatrixVector(A, v)
        lda = np.inner(w, v)
        k+=1
        delta_error = w - lda*v
        if np.linalg.norm(delta_error, ord=2) < matrix_size * epsilon:
            print("Iteration: {0}".format(k))
            return True, lda, v
    return False

def main():
    """ Generate my own rare matrix """
    #A = generateRareMatrix(600, 0.7)
    #matrix_size = len(A)
    #writeToFileGeneratedMatrix(A, 'D:\work\Anul3\SEM2\CN\GIT\Calcul_Numeric\Tema6\m_rar_sim_2018_2.txt')

    """ Generate the rare structures from the matrix in the input file """
    matrix_size, A = readFromFile('D:\work\Anul3\SEM2\CN\GIT\Calcul_Numeric\Tema6\m_rar_sim_2018.txt')

    """ Check the symmetry of the matrix read from file """
    if(checkMatrixSymmetry(A)):
        print("The matrix read from file is symmetric!\n")
        """ Determinarea celei mai mari valori proprii a matricei A si a unui vector propriu asociat """
        result = get_eigenvalues(matrix_size, A)
        if (result[0] == True):
            print("Valoarea proprie: {0}".format(result[1]))
            print("Vectorul asociat: {0}".format(result[2]))
    else:
        print("The matrix is not symmetric!\n")

    """ Case 3: p > n """
    """ A = readSVD("D:\work\Anul3\SEM2\CN\GIT\Calcul_Numeric\Tema6\SVD.txt")
    U, s, V = np.linalg.svd(A)
    print('Valorile singulare ale matricii: {}'.format([val for val in s]))
    print('Numarul de valori singulare din descompunere (rangul): {}'.format(len(s)))
    print('Rangul matricii: {}'.format(np.linalg.matrix_rank(A)))
    singular_vals = [val for val in s if val > 0]
    print('Numarul de conditionare al matricii: {}'.format(max(singular_vals) / min(singular_vals)))
    print('Numarul de conditionare al matricii folosind numpy: {0}'.format(np.linalg.cond(A)))
    B = np.linalg.pinv(A)
    print("Pseudoinversa Moore-Penrose a matricei A: \n{0}".format(B))
    b = np.random.random(5)
    x = np.dot(B, b)
    # The solution is the Moore-Penrose pseudoinverse, multiplied by b 
    print("Solution to Ax = b: {0}".format(x))
    s_val = int(input('s = '))
    A_s = np.zeros(A.shape)
    for i in range(s_val):
        col_u = U[:, i]
        col_v = V[:, i]
        A_s += s[i] * np.dot(col_u.reshape(len(col_u), 1), col_v.reshape(1, len(col_v)))
    print("Value of matrix A_s: \n")
    print(A_s)
    print('Norma ceruta: {}'.format(np.linalg.norm(A - A_s, np.inf))) """

if __name__ == '__main__':
    main()

