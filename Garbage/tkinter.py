from Tkinter import Tk, Label, Button, Entry, IntVar, END, W, E, Canvas, CENTER
import numpy as np
import tkMessageBox
import sys

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
    if boolean_value == False:
        return 0
    return 1

def check_multiplication_associativity():
    y = z = machine_precision_sum()[0]
    x = 1.0
    boolean_value = ((x * y) * z == x * (y * z))
    while boolean_value:
        x = x * 1.1
        boolean_value = ((x * y) * z == x * (y * z))
    return x, (x * y) * z, x * (y * z)

""" Problema 3 """
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

class Processor:
    def __init__(self, master):
        self.master = master
        master.geometry('500x500')
        master.title("Operation Processor")

        """ Prima problema """
        self.prob1 = Label(master, text="PROBLEMA 1:", font=("Helvetica bold", 10))
        self.prob1.pack()
        self.add_button = Button(master, text="Determine smallest number", fg = "blue", bg="white", command=lambda: self.update("prec"))
        self.min_no = IntVar()
        self.min_label = Label(master, textvariable=self.min_no)
        self.min_label.pack()

        """ A doua problema """
        self.prob2 = Label(master, text="PROBLEMA 2.1:", font=("Helvetica bold", 10))
        self.prob2.pack()
        self.add_button2 = Button(master, text="Check sum associativity", fg="blue", bg="white", command=lambda: self.update("assoc"))
        self.min_no2 = IntVar()
        self.min_label2 = Label(master, textvariable=self.min_no2)
        self.min_label2.pack()

        """ Al doilea subpunct - Problema 2"""
        self.prob22 = Label(master, text="PROBLEMA 2.2:", font=("Helvetica bold", 10))
        self.prob22.pack()
        self.add_button22 = Button(master, text="Determine element", fg="blue", bg="white", command=lambda: self.update("el"))
        self.min_no3 = IntVar()
        self.min_label3 = Label(master, textvariable=self.min_no3)
        self.min_label3.pack()

        """ A treia problema """
        self.prob3 = Label(master, text="PROBLEMA 3:", font=("Helvetica bold", 10))
        self.prob3.pack()
        self.prob31 = Label(master, text="Input n value:", font=("Helvetica bold", 10))
        self.prob31.pack()
        self.n_min = Entry(master)
        self.n_min.pack()

        # LAYOUT
        self.prob1.grid(row=0, column=0, sticky=W, pady = (20, 20))
        self.add_button.grid(row=0, column=1)
        self.min_label.grid(row=0, column=2)

        self.prob2.grid(row=2, column=0, sticky=W, pady = (20, 20))
        self.add_button2.grid(row=2, column=1)
        self.min_label2.grid(row=2, column=2)

        self.prob22.grid(row=4, column=0, sticky=W, pady = (20, 20))
        self.add_button22.grid(row=4, column=1)
        self.min_label3.grid(row=4, column=2)

        self.prob3.grid(row=6, column=0, sticky=W, pady = (5, 5))
        self.prob31.grid(row=7, column=0, sticky=W)
        self.n_min.grid(row=7, column=1, sticky=W)

    def update(self, method):
        if method == "prec":
            self.min_no.set(machine_precision_sum()[1])
        if method == "assoc":
            self.min_no2.set(check_non_associativity())
        if method == "el":
            self.min_no3.set(check_multiplication_associativity()[0])

if __name__ == "__main__":
    root = Tk()
    my_gui = Processor(root)
    root.mainloop()