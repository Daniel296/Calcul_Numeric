import numpy as np

def add_a_diagonal(A, a):
    for index in range(0, len(A)):
        A[index][index] += a
    return A;

def compute_V_1_1(V_0, B):
    """ Metoda Schultz de constructie a sirului / Algoritmul Hotelling - Bodewig
        sau metoda iterativa a hiper-puterii """
    C = np.dot(B, V_0)
    C = add_a_diagonal(C, 2)
    V_1 = np.dot(V_0, C)
    return V_1;

def compute_V_1_2(V_0, B):
    """ Metoda 2 propusa de Li si Li """
    C = np.dot(B, V_0)
    C = add_a_diagonal(C, 3)
    AV = np.dot(B, V_0)
    AV = np.dot(AV, C)
    AV = add_a_diagonal(AV, 3)
    V_1 = np.dot(V_0, AV)

    return V_1;

def compute_V_1_3(V_0, B):
    """ Metoda 3 propusa de Li si Li """
    C = np.dot(V_0, B)
    C = add_a_diagonal(C, 3)
    C = np.dot(C, C)
    P_1 = np.dot(V_0, B)
    P_1 = add_a_diagonal(P_1, 1)
    P_2 = np.dot(P_1, C)
    P_2 = np.divide(P_2, 4)
    V_1 = add_a_diagonal(P_2, 1)
    V_1 = np.dot(V_1, V_0)

    return V_1;

def determine_inverse(n, epsilon, k_max, A, index_metoda):
    norma_product = np.linalg.norm(A, np.inf) * np.linalg.norm(A, 1)
    A_transposed = np.transpose(A)
    initial_matrix = np.divide(A_transposed, norma_product)
    V_0 = initial_matrix
    V_1 = initial_matrix
    k = 0
    V_limit = pow(10, 10)
    """ Calculeaza o singura data in program B = (-A) """
    B = np.dot(A, -1)
    while True:
        V_0 = V_1
        if index_metoda == 1:
            V_1 = compute_V_1_1(V_0, B)
        if index_metoda == 2:
            V_1 = compute_V_1_2(V_0, B)
        if index_metoda == 3:
            V_1 = compute_V_1_3(V_0, B)
        if index_metoda not in [1, 2, 3]:
            print("Metoda neimplementata!\n")
            return False, None, None
        delta_V = np.linalg.norm(np.subtract(V_1, V_0))
        k += 1
        if (delta_V < epsilon or k > k_max or delta_V > V_limit):
            break
    if delta_V < epsilon:
        return True, k, V_1
    if k > k_max:
        print("Nu s-a reusit aproximarea inversei. Numar de iteratii depasit!\n")
        return False, k, None
    if delta_V > V_limit:
        print("Nu s-a reusit aproximarea inversei. Sirul construit nu converge!\n")
        return False, k, None

def induction_by_n(n, epsilon, k_max, A, index_metoda):
    raport = 0
    old_raport = 0
    if n < 2:
        print("Value of n too small for induction!\n")
        return False, None
    """ Determine the aferent smaller matrix starting from A """
    temp_A = np.zeros(shape=(2, 2))
    for line in range(0, 2):
        for col in range(0, 2):
            temp_A[line][col] = A[line][col]
    """ Determine the inverse of that temporary matrix """
    found, k, inverse_temp_A = determine_inverse(2, epsilon, k_max, temp_A, index_metoda)
    if found == False:
        print("Couldn't determine inverse of small matrix!\n")
        return False, None
    """ If the size is only 2, return it, else start the induction """
    if n == 2:
        return True, inverse_temp_A
    else:
        for size in range(3, n):
            old_raport = raport
            """ Save the last inverse """
            inverse_old_matrix = inverse_temp_A
            """ Construct the matrix for the next size """
            temp_A = np.zeros(shape=(size, size))
            for line in range(0, size):
                for col in range(0, size):
                    temp_A[line][col] = A[line][col]
            """ Determine the inverse for the new matrix """
            found, k, inverse_temp_A = determine_inverse(size, epsilon, k_max, temp_A, index_metoda)
            if found == False:
                print("Couldn't determine inverse of small matrix!\n")
                return False, None
            for line in range(0, len(inverse_old_matrix)):
                for column in range(0, len(inverse_old_matrix)):
                    if inverse_temp_A[line][column] - inverse_old_matrix[line][column] > epsilon:
                        print("Didn't find a pattern. Matrices are not alike!\n")
                        return False, None
                    for index in range(0, len(inverse_old_matrix)):
                        if inverse_temp_A[index][-1] != 0 and inverse_temp_A[index][-2] != 0:
                            raport = inverse_temp_A[index][-1] / inverse_temp_A[index][-2]
                            if raport != old_raport and old_raport != 0:
                                print("Didn't find a pattern. Ratio is not the same!\n")
                                return False, None
        inverse = np.zeros(shape=(n, n))
        """ First n-1 elements on line and column will be equal to the (n-1, n-1) last obtained inverse """
        for line in range(0, len(inverse_temp_A)):
            for col in range(0, len(inverse_temp_A)):
                inverse[line][col] += inverse_temp_A[line][col]
        """ Deduce the elements on the last column based on the ratio that is common """
        for line in range(0, n - 1):
            inverse[line][-1] =  inverse[line][-2] * raport
        inverse[n-1][n-1] = inverse[n-2][-2]

    return True, inverse

def readFromFile(file_path):
    fd = open(file_path, 'rb')
    """ Get matrix size """
    matrix_size = int(fd.readline().strip())
    """ Skip the white line """
    fd.readline()
    """ Get epsilon """
    epsilon_power = int(fd.readline().strip())
    epsilon = pow(10, -epsilon_power)
    """ Skip the white line """
    fd.readline()
    """ Get k_max """
    k_max = int(fd.readline().strip())
    """ Skip the white line """
    fd.readline()
    """ Initialize matrix A """
    A = np.zeros(shape=(matrix_size, matrix_size))
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
            col_index+=1
        line_index+=1

    return matrix_size, epsilon, k_max, A

def readFromFileAdapted(file_path):
    fd = open(file_path, 'rb')
    """ Get number of lines """
    lines = int(fd.readline().strip())
    """ Skip the white line """
    fd.readline()
    """ Get number of columns """
    cols = int(fd.readline().strip())
    """ Skip the white line """
    fd.readline()
    """ Get epsilon """
    epsilon_power = int(fd.readline().strip())
    epsilon = pow(10, -epsilon_power)
    """ Skip the white line """
    fd.readline()
    """ Get k_max """
    k_max = int(fd.readline().strip())
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
            col_index+=1
        line_index+=1

    return lines, cols, epsilon, k_max, A

def main():
    n, epsilon, k_max, A = readFromFile("input.txt")
    index_metoda = 3
    found, k, inverse = determine_inverse(n, epsilon, k_max, A, index_metoda)
    np.set_printoptions(suppress=True)
    if found == True:
        #print("Convergenta la iteratia:{0}\n".format(k))
        #print(inverse)
        """ Compute norm """
        mul = np.dot(A, inverse)
        diff = add_a_diagonal(mul, -1)
        norma = np.linalg.norm(diff)
        #print("\n Norma: {0}\n".format(norma))
    """ Cerinta 3 """
    """ n, epsilon, k_max, A = readFromFile("input1.txt")
    found, inverse = induction_by_n(n, epsilon, k_max, A, index_metoda)
    if found == True:
        print(inverse)
        found, k, inverse = determine_inverse(n, epsilon, k_max, A, index_metoda)
        np.set_printoptions(suppress=True)
        print("Inverse by function: \n")
        print(inverse) """
    """ Bonus """
    no_lines, no_cols, epsilon, k_max, A = readFromFileAdapted("input2.txt")
    if np.linalg.matrix_rank(A) == no_lines:
        print(""" Determine the right inverse!\n """)
        found, k, inverse = determine_inverse(n, epsilon, k_max, A, 2)
        if found == True:
            print(inverse)
    if np.linalg.matrix_rank(A) == no_cols:
        print(""" Determine the left inverse!\n """)
        found, k, inverse = determine_inverse(n, epsilon, k_max, A, 2)
        if found == True:
            print(inverse)
if __name__ == '__main__':
    main()
