#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: zlel_main.py
    :synopsis:

.. moduleauthor:: Aitor Lavin (aitorlavin02@gmail.com)
                  Unax Arregi (arregiunax@gmail.com)


"""

import numpy as np
import sys
import matplotlib.pyplot as plt

if __name__ == "zlel.zlel_p2":
    import zlel.zlel_p1 as zl1
else:
    import zlel_p1 as zl1

"IRAULIA"


def get_iraulia(incidence_matrix):
    """
    np array bat ematen du matrizearen irauliarena.

    Parameters
    ----------
    incidence_matrix : np array-a zirkuituaren intzidentzia matrizearena.

    Returns
    -------
    np_iraulia : np array-a matrize irauliarena.

    """

    num_filas = len(incidence_matrix)
    num_columnas = len(incidence_matrix[0])
    iraulia = [[0] * num_filas for _ in range(num_columnas)]
    for i in range(num_filas):
        for j in range(num_columnas):
            iraulia[j][i] = incidence_matrix[i][j]
    np_iraulia = np.array(iraulia)
    return (np_iraulia)


"M N ETA Us LORTZEKO BEHAR DE FUNTZIOA"


def getPosition(x, cir_el_extended):
    """
    cir_ctr_extended-ko x balioko elementua cir_el_extendedeko
    i-koaren berdina bada i posizioa bueltatzen du.

    Parameters
    ----------
    x : String bat cir_ctr_extended-eko elementu batena.
    cir_el_extended : cir_el array luzatua.

    Returns
    -------
    i : integer-a cir_el_extendededeko elementu baten posizioa
    adierazten duena.

    """

    for i in range(0, np.size(cir_el_extended)):
        if cir_el_extended[i].lower() == x.lower():
            return i


"M, N ETA Us LORTZEKO FUNTZIOA"


def getMatrixes(b, cir_el_extended, cir_val_extended, cir_ctr_extended):
    """
    M, N ETA Us matrizeak lortzeko funtzioa

    Parameters
    ----------
    b : integer
        Gure zirkuituko adar kopurua.
    cir_el_extended : cir_el array luzatua.
    cir_val_extended : cir_val array luzatua.
    cir_ctr_extended : cir_ctr array luzatua.

    Returns
    -------
    lista : Zerrenda bat M, N eta Us matrizeak dituena.

    """

    M = np.zeros((b, b), dtype=float)
    N = np.zeros((b, b), dtype=float)
    Us = np.zeros((b, 1), dtype=float)
    i = 0
    while i < b:

        if cir_el_extended[i][0].lower() == "r":
            M[i][i] = 1
            N[i][i] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "v":
            M[i][i] = 1
            Us[i] = float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "i":
            N[i][i] = 1
            Us[i] = float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "a":
            if "in" in cir_el_extended[i]:
                M[i][i] = 1
                N[i+1][i] = 1
        if cir_el_extended[i][0].lower() == "e":
            j = np.where(np.char.lower(cir_el_extended) == np.char.lower(
                cir_ctr_extended[i]))[0]
            M[i][i] = 1
            M[i][j] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "g":
            j = np.where(np.char.lower(cir_el_extended) == np.char.lower(
                cir_ctr_extended[i]))[0]
            N[i][i] = 1
            M[i][j] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "f":
            j = getPosition(cir_ctr_extended[i], cir_el_extended)
            N[i][i] = 1
            N[i][j] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "h":
            j = getPosition(cir_ctr_extended[i], cir_el_extended)
            M[i][i] = 1
            N[i][j] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "b":
            M[i][i] = 1
            Us[i] = float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "y":
            N[i][i] = 1
            Us[i] = float(cir_val_extended[i][0])
        i += 1
    lista = []
    lista.append(M)
    lista.append(N)
    lista.append(Us)
    return lista


"w(t) MATRIZEA"


def w_matrizea(nodes, b):
    """
    w(t) MATRIZEA lortzeko funtzioa.

    Returns
    -------
    matriz_np : np array bat w(t) matrizearena.

    """

    elementos = []
    for i in nodes[1:]:
        elementos.append(["e" + str(i)])

    for i in range(b-1):
        elementos.append(["v" + str(i+1)])
    elementos.append(["v" + str(b)])

    for i in range(b):
        elementos.append(["i" + str(i+1)])

    matriz_np = np.array(elementos)
    return (matriz_np)


"IDENTITATE MATRIZEA"


def get_identitate_matrizea(M, iraulia):
    """
    Matrize irauliaren errenkada kopurua eta M matrizearen zutabe kopurua duen
    matrizea bueltatzen du.

    Parameters
    ----------
    M : np array-a M matrizearena.
    iraulia : np array-a matrize irauliarena.

    Returns
    -------
    Identitate matrizea.

    """

    rows_iraulia, _ = iraulia.shape
    _, cols_M = M.shape
    return np.eye(rows_iraulia, cols_M)


"T MATRIZEA"


def get_T(incidence_matrix, iraulia, M, N, Us, b, n):
    """
    T matrizea bueltatzen du.

    Parameters
    ----------
    incidence_matrix : Intzidentzia matrizea
    iraulia : np array-a matrize irauliarena.
    M : np array-a M matrizearena.
    N : np array-a N matrizearena.
    Us : np array-a Us matrizearena.
    b : integer-a adar kopuruarena.
    n : Integer-a nodo kopuruarena.

    Returns
    -------
    T : np array-a T matrizearena.

    """
    T = np.zeros((n+2*b-1, n+2*b-1), dtype=float)
    T[:n-1, n-1+b:] = incidence_matrix
    T[n-1:n-1+b, :n-1] = -iraulia
    T[n-1:n-1+b, n-1:n-1+b] = np.eye(b)
    T[n-1+b:, n-1:n-1+b] = M
    T[n-1+b:, n-1+b:] = N

    return T


"U MATRIZEA"


def get_U(incidence_matrix, b, n, Us):
    """
    Funtzio honek U matrizea lortzeko balio du.

    Parameters
    ----------
    incidence_matrix : Intzidentzia matrizea
    b : integer-a adar kopuruarena.
    n : Integer-a nodo kopuruarena.
    Us : np array-a Us matrizearena.

    Returns
    -------
    U : np array-a U matrizearena.

    """

    U = np.zeros((n-1+2*b, 1), dtype=float)
    U[n-1+b:n-1+2*b] = Us

    return U


"SINGULAR MATRIX ERROREA"


def SingularMatrix(T, U):
    """
    Singular matrix errorea ematen bada,hau da, T-ren determinantea 0 ez bada,
    errorea gertatu dela adierazten du.

    Parameters
    ----------
    T : np array-a T matrizearena.
    U : np array-a U matrizearena.

    Returns
    -------
    None.

    """

    if np.linalg.det(T) == 0:
        sys.exit("Error solving Tableau equations, check if det(T) != 0.")


def print_solution(sol, b, n):
    """ This function prints the solution with format.

        Args:
            | sol: np array with the solution of the Tableau equations
            | (e_1,..,e_n-1,v_1,..,v_b,i_1,..i_b)
            | b: # of branches
            | n: # of nodes

    """

    # The instructor solution needs to be a numpy array of numpy arrays of
    # float. If it is not, convert it to this format.
    if sol is None:
        pass
    else:
        if sol.dtype == np.float64:
            np.set_printoptions(sign=' ')  # Only from numpy 1.14
            tmp = np.zeros([np.size(sol), 1], dtype=float)
            for ind in range(np.size(sol)):
                tmp[ind] = np.array(sol[ind])
            sol = tmp
        print("\n========== Nodes voltage to reference ========")
        for i in range(1, n):
            print("e" + str(i) + " = ", "[{:10.9f}]".format(sol[i-1][0]))
        print("\n========== Branches voltage difference ========")
        for i in range(1, b+1):
            print("v" + str(i) + " = ", "[{:10.9f}]".format(sol[i+n-2][0]))
        print("\n=============== Branches currents ==============")
        for i in range(1, b+1):
            print("i" + str(i) + " = ", "[{:10.9f}]".format(sol[i+b+n-2][0]))

        print("\n================= End solution =================\n")


def build_csv_header(tvi, b, n):
    """ This function build the csv header for the output files.
        First column will be v or i if .dc analysis or t if .tr and it will
        be given by argument tvi.
        The header will be this form,
        t/v/i,e_1,..,e_n-1,v_1,..,v_b,i_1,..i_b

    Args:
        | tvi: "v" or "i" if .dc analysis or "t" if .tran
        | b: # of branches
        | n: # of nodes

    Returns:
        header: The header in csv format as string
    """
    header = tvi
    for i in range(1, n):
        header += ",e" + str(i)
    for i in range(1, b+1):
        header += ",v" + str(i)
    for i in range(1, b+1):
        header += ",i" + str(i)
    return header


def save_as_csv(b, n, filename, cir_el_extended, cir_val_extended,
                cir_ctr_extended, A, At, cir_nd, cir_nd_extended,
                incidence_matrix, nodes, a, hasiera, amaiera, pausua,
                sorgailua):
    """ This function generates a csv file with the name filename.
        First it will save a header and then, it loops and save a line in
        csv format into the file.

    Args:
        | b: # of branches
        | n: # of nodes
        | filename: string with the filename (incluiding the path)
    """
    # Sup .tr
    header = build_csv_header("t", b, n)
    with open(filename, 'w') as file:
        print(header, file=file)
        # Get the indices of the elements corresponding to the sources.
        # The freq parameter cannot be 0 this is why we choose cir_tr[0].
        t = hasiera
        while t < amaiera:
            if a == ".tr":
                m = op_tr(b, cir_el_extended, cir_val_extended,
                          cir_ctr_extended, t)
            elif a == ".dc":
                m = op_dc(b, cir_el_extended, cir_val_extended,
                          cir_ctr_extended, t, sorgailua)
            M = m[0]
            N = m[1]
            Us = m[2]
            T = get_T(A, At, M, N, Us, b, n)
            U = get_U(A, b, n, Us)
            SingularMatrix(T, U)
            zl1.erroreak(cir_nd_extended, cir_el_extended, cir_ctr_extended,
                         cir_val_extended, incidence_matrix, nodes, b)
            sol = np.linalg.solve(T, U)
            # for t in tr["start"],tr["end"],tr["step"]
            # Recalculate the Us for the sinusoidal sources

            # sol = np.full(2*b+(n-1), t+1, dtype=float)
            # Inserte the time
            sol = np.insert(sol, 0, t)
            # sol to csv
            sol_csv = ','.join(['%.9f' % num for num in sol])
            print(sol_csv, file=file)
            t = t + pausua


def plot_from_cvs(filename, x, y, title):
    """ This function plots the values corresponding to the x string of the
        file filename in the x-axis and the ones corresponding to the y
        string in the y-axis.
        The x and y strings must mach with some value of the header in the
        csv file filename.

    Args:
        | filename: string with the name of the file (including the path).
        | x: string with some value of the header of the file.
        | y: string with some value of the header of the file.

    """
    data = np.genfromtxt(filename, delimiter=',', skip_header=0,
                         skip_footer=1, names=True)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(data[x], data[y], color='r', label=title)
    ax1.set_xlabel(x)
    ax1.set_ylabel(y)
    plt.show()


def op_dc(b, cir_el_extended, cir_val_extended, cir_ctr_extended, t,
          sorgailua):
    M = np.zeros((b, b), dtype=float)
    N = np.zeros((b, b), dtype=float)
    Us = np.zeros((b, 1), dtype=float)
    i = 0
    while i < b:

        if cir_el_extended[i][0].lower() == "r":
            M[i][i] = 1
            N[i][i] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "v":
            M[i][i] = 1
            Us[i] = float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "i":
            N[i][i] = 1
            Us[i] = float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "a":
            if "in" in cir_el_extended[i]:
                M[i][i] = 1
                N[i+1][i] = 1
        if cir_el_extended[i][0].lower() == "e":
            j = np.where(np.char.lower(cir_el_extended) == np.char.lower(
                cir_ctr_extended[i]))[0]
            M[i][i] = 1
            M[i][j] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "g":
            j = np.where(np.char.lower(cir_el_extended) == np.char.lower(
                cir_ctr_extended[i]))[0]
            N[i][i] = 1
            M[i][j] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "f":
            j = getPosition(cir_ctr_extended[i], cir_el_extended)
            N[i][i] = 1
            N[i][j] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "h":
            j = getPosition(cir_ctr_extended[i], cir_el_extended)
            M[i][i] = 1
            N[i][j] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "b":
            M[i][i] = 1
            Us[i] = float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "y":
            N[i][i] = 1
            Us[i] = float(cir_val_extended[i][0])
        i += 1
    j = np.where(np.char.lower(cir_el_extended) == np.char.lower(sorgailua))[0]
    Us[j] = t
    lista = []
    lista.append(M)
    lista.append(N)
    lista.append(Us)
    return lista


def op_tr(b, cir_el_extended, cir_val_extended, cir_ctr_extended, t):
    M = np.zeros((b, b), dtype=float)
    N = np.zeros((b, b), dtype=float)
    Us = np.zeros((b, 1), dtype=float)
    i = 0
    while i < b:

        if cir_el_extended[i][0].lower() == "r":
            M[i][i] = 1
            N[i][i] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "v":
            M[i][i] = 1
            Us[i] = float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "i":
            N[i][i] = 1
            Us[i] = float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "a":
            if "in" in cir_el_extended[i]:
                M[i][i] = 1
                N[i+1][i] = 1
        if cir_el_extended[i][0].lower() == "e":
            j = np.where(np.char.lower(cir_el_extended) == np.char.lower(
                cir_ctr_extended[i]))[0]
            M[i][i] = 1
            M[i][j] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "g":
            j = np.where(np.char.lower(cir_el_extended) == np.char.lower(
                cir_ctr_extended[i]))[0]
            N[i][i] = 1
            M[i][j] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "f":
            j = getPosition(cir_ctr_extended[i], cir_el_extended)
            N[i][i] = 1
            N[i][j] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "h":
            j = getPosition(cir_ctr_extended[i], cir_el_extended)
            M[i][i] = 1
            N[i][j] = -float(cir_val_extended[i][0])
        if cir_el_extended[i][0].lower() == "b":
            anp = float(cir_val_extended[i][0])
            maiz = float(cir_val_extended[i][1])
            fas = float(cir_val_extended[i][2])  # Conversión a float
            m = 2 * np.pi * maiz
            f = (np.pi / 180) * fas  # Ahora sí funciona
            M[i, i] = 1
            Us[i] = anp * np.sin(m * t + f)
        if cir_el_extended[i][0].lower() == "y":
            anp = cir_val_extended[i][0]
            maiz = float(cir_val_extended[i][1])
            fas = float(cir_val_extended[i][2])
            m = 2*np.pi*maiz
            f = float((np.pi/180)*fas)
            N[i, i] = 1
            Us[i] = anp * np.sin(m*t+f)
        i += 1
    lista = []
    lista.append(M)
    lista.append(N)
    lista.append(Us)
    return lista


def tr_dc_parametroak(filename):
    """ Esta función carga los valores de un archivo .cir
    y extrae los parámetros.

    Args:
        filename: ruta del archivo de circuito.

    Returns:
        A tuple with (sorgailua, hasiera, amaiera, pausua).
    """
    sorgailua = None
    hasiera, amaiera, pausua = None, None, None
    a = None  # Tipo de análisis (.tr o .dc)
    pr = None  # Printa egin behar den adierazi
    op= None  # Operazio puntua egin behar den adierazi
    # Leer el archivo .cir
    with open(filename, 'r') as file:
        lines = file.readlines()

        for line in lines:
            tokens = line.strip().split()

            if len(tokens) < 1:
                continue  # Omitir líneas vacías

            # Detectar el análisis transitorio (.tr)
            if tokens[0].lower() == ".tr" or tokens[0].lower() == ".dc":
                a = tokens[0].lower()
                hasiera = float(tokens[5])  # Tiempo inicial
                amaiera = float(tokens[6])  # Tiempo final
                pausua = float(tokens[7])   # Paso de tiempo
                sorgailua = str(tokens[8])

            if tokens[0].lower() == ".pr":
                pr = ".pr"

            if tokens[0].lower() == ".op":
                op = ".op"

    return sorgailua, hasiera, amaiera, pausua, a, pr, op


"""
https://stackoverflow.com/questions/419163/what-does-if-name-main-do
https://stackoverflow.com/questions/19747371/
python-exit-commands-why-so-many-and-when-should-each-be-used
"""


if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "cirs/all/1_zlel_serial_YI_VI.cir"

    "1 PRAKTIKAKO BALIOAK LORTU"

    cir_el, cir_nd, cir_val, cir_ctr = zl1.cir_parser(filename)
    modified_results = zl1.modify_function(cir_el, cir_nd, cir_val, cir_ctr)
    cir_el_extended = np.array(modified_results[0])
    cir_nd_extended = np.array(modified_results[1])
    cir_ctr_extended = np.array(modified_results[3])
    cir_val_extended = np.array(modified_results[2])
    nodes = zl1.get_nodes(cir_nd)
    b = zl1.get_branches(cir_el)
    n = zl1.get_number_of_nodes(nodes)
    el_num = zl1.get_elements(cir_el)
    incidence_matrix = zl1.get_incidence_matrix(n, b, nodes, cir_nd_extended)
    iraulia = get_iraulia(incidence_matrix)

    "INTZIDENTZI MATRIZE MURRIZTUA"

    A = zl1.get_incidence_murriztua(incidence_matrix, n)

    "INTZIDENTZI MATRIZE MURRIZTUAREN IRAULIA"

    At = get_iraulia(A)

    "PRAKTIKA HONETAKO FUNTZIOEN BALIOAK LORTU"

    m = getMatrixes(b, cir_el_extended, cir_val_extended, cir_ctr_extended)
    M = m[0]
    N = m[1]
    Us = m[2]

    "w(T)"

    w = w_matrizea(nodes, b)

    "IDENTITATE MATRIZEA"

    identitate_matrizea = get_identitate_matrizea(M, At)

    "T ETA U"

    T = get_T(A, At, M, N, Us, b, n)
    U = get_U(A, b, n, Us)

    "ERROREAK AURKITU"

    zl1.erroreak(cir_nd_extended, cir_el_extended, cir_ctr_extended,
                 cir_val_extended, incidence_matrix, nodes, b)

    "PRAKTIKA HONETAKO ERROREA (SINGULAR MATRIX)"

    SingularMatrix(T, U)

    "SOLUZIOA"

    sol = np.linalg.solve(T, U)
    
    "sIMULAZIOAK"
    sorgailua, hasiera, amaiera, pausua, a, pr = tr_dc_parametroak(filename)
    # Definir archivo de salida

    if pr == ".pr":
        print_solution(sol, b, n)

    if a == ".tr" or a == ".dc":
        # Ejecutar simulación y guardar resultados en CSV
        output_csv = "resultado.csv"
        save_as_csv(
            b, n, output_csv, cir_el_extended, cir_val_extended,
            cir_ctr_extended, A, At, cir_nd, cir_nd_extended,
            incidence_matrix, nodes, a, hasiera, amaiera, pausua, sorgailua
            )

        # Graficar resultado de la simulación
        plot_from_cvs(output_csv, "t", "e1", "Voltaje en nodo 1")
