
def readFromFile(file_path):
    fd = open(file_path, 'rb')
    values = {}
    """ Read the values of the function """
    while True:
        line = fd.readline()
        array = line.split()
        if len(array) < 2:
            break
        values[float(array[0])] = float(array[1])
    """ Get value of n """
    line = fd.readline().strip()
    n_array = line.split("=")
    if len(n_array) < 2:
        return False, None, None, None, None, None
    n = int(n_array[1])
    """ Get value of x0 """
    line = fd.readline().strip()
    x0_array = line.split("=")
    if len(x0_array) < 2:
        return False, None, None, None, None, None
    x0 = float(x0_array[1])
    """ Get value of xn """
    line = fd.readline().strip()
    xn_array = line.split("=")
    if len(xn_array) < 2:
        return False, None, None, None, None, None
    xn = float(xn_array[1])
    """ Get value of x """
    line = fd.readline().strip()
    x_array = line.split("=")
    if len(x_array) < 2:
        return False, None, None, None, None, None
    x = float(x_array[1])

    return True, values, n, x0, xn, x

def progresiveNewton(values, n, x0, xn, x):
    if xn < x0:
        return False, None
    h = (xn - x0)/(n-1)
    """ Determine the interpolation points """
    x_interpol = [0] * (n)
    x_interpol[0] = x0
    for i in range(1, n):
        x_interpol[i] = x0 + i * h
    """ Determine t """
    t = (x - x0)/h
    """ Start computing L """
    y = [0] * (n)
    y = [0] * (n)
    s = [0] * (n)
    y[0] = values[x0]
    s[0] = 1
    s[1] = t
    L = 0

    """ First step """
    for counter in range(1, n):
        y[counter] = values[values.keys()[counter + 1]] - values[values.keys()[counter]]

    """ Following steps """
    for pas in range(2, n):
        s[pas] = s[pas - 1] * (t - pas + 1) / pas

        z = [0] * (n-pas)
        for counter in range(0, n-pas):
            z[counter] = y[pas + counter] - y[pas + counter - 1]
        for counter in range(0, n-pas):
            y[pas + counter] = z[counter]

    """ Compute L """
    for counter in range(0, n):
        L += s[counter] * y[counter]

    return True, L

def main():
    check, values, n, x0, xn, x = readFromFile("input.txt")
    if check == True:
        check, result = progresiveNewton(values, n, x0, xn, x)
    if check == True:
        print(result)
    else:
        print("Couldn't find result!")

if __name__ == '__main__':
    main()
