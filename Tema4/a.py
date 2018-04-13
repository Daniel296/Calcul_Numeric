import numpy


def citeste_matrice(fisier):
    f = open(fisier)
    
    txt = f.readlines()
    n = int(txt[0])
    b = []
    matrice = []

    for i in range(n):
        matrice.append([])

    for i in range(2, n + 2):
        b.append(float(txt[i]))

    for index in range(n + 3, len(txt)):
        parsat = txt[index].replace(' ', '').split(',')
        elem = float(parsat[0])
        i = int(parsat[1])
        j = int(parsat[2])

        exista_deja = 0
        for coloana in range(0,len(matrice[i])):
            if matrice[i][coloana][1] == j:
                matrice[i][coloana][0] += elem
                exista_deja = 1

        if exista_deja == 0:
            if i == j:
                matrice[i].append([elem, j])
            else:
                matrice[i] = [[elem,j]] + matrice[i]
    return (n, b, matrice)


def verificare_diag(nr, matrice, n):
    for i in range(n):
        exista_diag = 0
        for j in range(len(matrice[i])):
            if i == matrice[i][j][1]:
                exista_diag = 1
                if abs(matrice[i][j][0]) <= epsilon:
                    print(nr, " Matricea nebuna!")
                    return False
        if exista_diag == 0:
            print(nr, " Matricea nebuna!")
            return False

    print(nr, " Matrice buna!")
    return True

def calculeaza(nr, matrice, b, n):
    #xgs = numpy.zeros(n)
    xgs = [0] * n

    for k in range(kmax):
        # print("k=",k)
        convergent = 1
        divergent = 1
        for i in range(n):

            # calculam noul element
            suma1 = suma2 = 0
            contor1 = 0
            while contor1 < len(matrice[i]) and  matrice[i][contor1][1] < i:
                suma1 += matrice[i][contor1][0] * xgs[matrice[i][contor1][1]]
                contor1 += 1

            contor1 += 1

            while contor1 < len(matrice[i]) and matrice[i][contor1][1] < n:
                suma2 += matrice[i][contor1][0] * xgs[matrice[i][contor1][1]]
                contor1 += 1

            ind_diag = 0
            while matrice[i][ind_diag][1] != i:
                ind_diag += 1
            elem = (b[i] - suma1 - suma2) / matrice[i][ind_diag][0]

            # verificam delta x
            # if abs(xgs[i] - elem) >= epsilon and abs(xgs[i] - elem) <= pow(10, 8):
            #     gata = 0

            if abs(xgs[i] - elem) >= epsilon:
                convergent = 0

            if abs(xgs[i] - elem) <= pow(10, 8):
                divergent = 0

            xgs[i] = elem

        #print("noul x:", xgs)
        if convergent or divergent:
            break


    print("k = ", k)
    if convergent and k != kmax - 1:
        # print("solutie: ", xgs)
        return xgs
    else:
        print("divergent!!")
        return None


def inmulteste_matrice_cu_vector(A, x, n):
    C = []

    for i in range(n):  #liniile lui A
        elem = 0
        for k in range(len(A[i])):
            index_x = A[i][k][1]
            elem += A[i][k][0] * x[index_x]
        if elem != 0:
            C.append(elem)

    return C



epsilon = pow(10, -5)
kmax = 10000

n1, b1, A1 = citeste_matrice("m_rar_2018_2.txt")
if verificare_diag(1, A1, n1):
    print("sorteaza matricea...")
    for i in range(n1):
        A1[i].sort(key=lambda x: x[1])
    print("sort done")
    
    xgs = calculeaza(1, A1, b1, n1)
    if xgs != None:
        C = inmulteste_matrice_cu_vector(A1, xgs, n1)
        #print("C: ")
        #print(C)

        for i in range(n1):
            C[i] -= b1[i]

        

        print("Norma: ")
        print(numpy.linalg.norm(C, numpy.inf))

        # print("diferenta:")
        # print(C)

        print("solutie: ", xgs)