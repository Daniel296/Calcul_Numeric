import sys
import numpy as np
from copy import deepcopy
from scipy.spatial import distance
from math import sqrt
from PyQt4 import QtCore, QtGui, uic

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
            results_copy[it] = results_copy[it + 1]
            results_copy[it + 1] = aux

        temp = (-1 * d[it] / cc[it])
        d[it + 1] = temp * d[it + 1] + e[it]
        if it < diag_size - 2:
            e[it + 1] = temp * e[it + 1] + f[it]
        if it < diag_size - 3:
            f[it] *= temp
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

""" UI Part """
qtCreatorFile = "D:\work\Anul3\SEM2\CN\GIT\Calcul_Numeric\Main_Frame2.ui"

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class MyApp(QtGui.QMainWindow, Ui_MainWindow):

    def __init__(self):
        super(MyApp, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.Result_2.clicked.connect(self.Problema1)
        self.ui.Result_tri.clicked.connect(self.Problema2)
        self.ui.V1.clicked.connect(self.Verificare1)
        self.ui.V2.clicked.connect(self.Verificare2)
        self.ui.V3.clicked.connect(self.Verificare3)

    def Problema1(self):
        self.ui.statement.setText("Sa se aproximeze solutia ecuatiei: Ax = b, folosindu-se Gauss cu pivotare partiala si metoda substitutiei inverse.\n")
        n = int(self.ui.n_val.toPlainText())
        path = self.ui.path.toPlainText()
        epsilon = self.ui.epsilon.toPlainText()
        epsilon = pow(10, (-1)*int(epsilon))
        fd = open(path, 'r')
        buf = fd.readlines()
        fd.close()
        A = np.zeros((n, n))
        for line in range(0, n):
            A[line] = np.array([float(i.strip()) for i in buf[line].split(' ')])
        b = np.array([float(i.strip()) for i in buf[-1].split(' ')])
        my_result = partial_gauss(n, epsilon, A, b)
        matrix = ''
        for i in range(0, len(my_result)):
            matrix += str(format(my_result[i], '.2f'))
            matrix += " "
        self.ui.Result.setText(matrix)

    def Problema2(self):
        self.ui.statement.setText("Sa se aproximeze solutia ecuatiei: Ax= b, folosindu-se Gauss cu pivotare partiala, adaptat pentru matrici tridiagonale.\n")
        n = int(self.ui.n_val.toPlainText())
        path = self.ui.path.toPlainText()
        epsilon = self.ui.epsilon.toPlainText()
        epsilon = pow(10, (-1) * int(epsilon))
        fd = open(path, 'r')
        buf = fd.readlines()
        fd.close()
        A = np.zeros((n, n))
        for line in range(0, n):
            A[line] = np.array([float(i.strip()) for i in buf[line].split(' ')])
        b = np.array([float(i.strip()) for i in buf[-1].split(' ')])
        Acopy_2 = deepcopy(A)
        bcopy_2 = deepcopy(b)
        a = [Acopy_2[i][i] for i in range(n)]
        b = [Acopy_2[i][i + 1] for i in range(n - 1)]
        c = [Acopy_2[i + 1][i] for i in range(n - 1)]
        d = deepcopy(bcopy_2)
        my_result = tri_diag_matrix_solver(a, b, c, d)
        matrix = ''
        for i in range(0, len(my_result)):
            matrix += str(format(my_result[i], '.2f'))
            matrix += " "
        self.ui.Result.setText(matrix)

    def Verificare1(self):
        self.ui.statement.setText("Sa se afiseze, spre verificare: || A_init * x_Gauss - b_init ||.\n")
        n = int(self.ui.n_val.toPlainText())
        path = self.ui.path.toPlainText()
        epsilon = self.ui.epsilon.toPlainText()
        epsilon = pow(10, (-1) * int(epsilon))
        fd = open(path, 'r')
        buf = fd.readlines()
        fd.close()
        A = np.zeros((n, n))
        for line in range(0, n):
            A[line] = np.array([float(i.strip()) for i in buf[line].split(' ')])
        b = np.array([float(i.strip()) for i in buf[-1].split(' ')])
        Acopy = deepcopy(A)
        bcopy = deepcopy(b)
        my_result = partial_gauss(n, epsilon, A, b)
        result_mul = matrixXvector(Acopy, my_result)
        result_sub = vectorMinusVector(result_mul, bcopy)
        norma = euclidianNorma(result_sub)
        self.ui.Result.setText(str(norma))

    def Verificare2(self):
        self.ui.statement.setText("Sa se afiseze, spre verificare: || x_Gauss - x_lib ||.\n")
        n = int(self.ui.n_val.toPlainText())
        path = self.ui.path.toPlainText()
        epsilon = self.ui.epsilon.toPlainText()
        epsilon = pow(10, (-1) * int(epsilon))
        fd = open(path, 'r')
        buf = fd.readlines()
        fd.close()
        A = np.zeros((n, n))
        for line in range(0, n):
            A[line] = np.array([float(i.strip()) for i in buf[line].split(' ')])
        b = np.array([float(i.strip()) for i in buf[-1].split(' ')])
        Acopy = deepcopy(A)
        bcopy = deepcopy(b)
        my_result = partial_gauss(n, epsilon, A, b)
        np_result = (np.linalg.solve(Acopy, bcopy))
        result_sub = vectorMinusVector(my_result, np_result)
        norma = euclidianNorma(result_sub)
        self.ui.Result.setText(str(norma))

    def Verificare3(self):
        self.ui.statement.setText("Sa se afiseze, spre verificare: || x_Gauss - A_lib_inv* b_init ||.\n")
        n = int(self.ui.n_val.toPlainText())
        path = self.ui.path.toPlainText()
        epsilon = self.ui.epsilon.toPlainText()
        epsilon = pow(10, (-1) * int(epsilon))
        fd = open(path, 'r')
        buf = fd.readlines()
        fd.close()
        A = np.zeros((n, n))
        for line in range(0, n):
            A[line] = np.array([float(i.strip()) for i in buf[line].split(' ')])
        b = np.array([float(i.strip()) for i in buf[-1].split(' ')])
        Acopy = deepcopy(A)
        bcopy = deepcopy(b)
        my_result = partial_gauss(n, epsilon, A, b)
        inverse = np.linalg.inv(Acopy)
        result_mul = matrixXvector(inverse, bcopy)
        result_sub = vectorMinusVector(my_result, result_mul)
        norma = euclidianNorma(result_sub)
        self.ui.Result.setText(str(norma))

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = MyApp()
    window.show()
    sys.exit(app.exec_())