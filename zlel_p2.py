#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: zlel_main.py
    :synopsis: zlel_p1 moduluan egindakoari zirkuitu erresistibo linealak
              (diodoa eta transistorea ez ditugu kontuan hartuko) ebazteko
              gaitasuna gehituko diogu (.OP). Halaber, programa .DC eta .TR
              (.TRAN spice-n) analisiak ere egiteko gai izango da.

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


def get_iraulia(incidence_matrix):
    """

    np array bat ematen du matrizearen irauliarena.

    Args
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


def getPosition(x, cir_el_extended):
    """

    cir_ctr_extended-ko x balioko elementua cir_el_extendedeko
    i-koaren berdina bada i posizioa bueltatzen du.

    Args
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


def getMatrixes(b, cir_el_extended, cir_val_extended, cir_ctr_extended):
    """

    M, N ETA Us matrizeak lortzeko funtzioa

    Args
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


def w_matrizea(nodes, b):
    """

    w(t) MATRIZEA lortzeko funtzioa.

    Args
    ----------
    nodes : np array-a zirkuituko nodoekin.
    b : Integer-a zirkuituko elementu kopuruarena.

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


def get_identitate_matrizea(M, iraulia):
    """

    Matrize irauliaren errenkada kopurua eta M matrizearen zutabe kopurua duen
    matrizea bueltatzen du.

    Args
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


def get_T(incidence_matrix, iraulia, M, N, Us, b, n):
    """

    T matrizea bueltatzen du.

    Args
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


def get_U(incidence_matrix, b, n, Us):
    """

    Funtzio honek U matrizea lortzeko balio du.

    Args
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


def SingularMatrix(T, U):
    """
    Singular matrix errorea ematen bada,hau da, T-ren determinantea 0 ez bada,
    errorea gertatu dela adierazten du.

    Args
    ----------
    T : np array-a T matrizearena.
    U : np array-a U matrizearena.

    Returns
    -------
    String bat errorea azalduz.

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
    if sol is None:
        pass
    else:
        if sol.dtype == np.float64:
            np.set_printoptions(sign=' ')
            tmp = np.zeros([np.size(sol), 1], dtype=float)
            for ind in range(np.size(sol)):
                tmp[ind] = np.array(sol[ind])
            sol = tmp
        print("\n========== Nodes voltage to reference ========")
        for i in range(1, n):
            value = float(sol[i-1][0])
            print("e" + str(i) + " = ",
                  "[{:10.9f}]".format(0.0 if abs(value) < 1e-12 else value))
        print("\n========== Branches voltage difference ========")
        for i in range(1, b+1):
            value = sol[i+n-2][0]
            value = 0.0 if abs(value) < 1e-12 else value
            print("v" + str(i) + " = ", "[{:10.9f}]".format(value))
        print("\n=============== Branches currents ==============")
        for i in range(1, b+1):
            value = sol[i+b+n-2][0]
            print("i" + str(i) + " = ",
                  "[{:10.9f}]".format(0.0 if abs(value) < 1e-12 else value))

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
    header = build_csv_header("t", b, n)
    with open(filename, 'w') as file:
        print(header, file=file)
        t = hasiera
        num_steps = int(round((amaiera - hasiera) / pausua)) + 1
        for i in range(num_steps):
            t = hasiera + i * pausua
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
            sol = np.insert(sol, 0, t)
            sol_csv = ','.join(['%.9f' % num for num in sol])
            print(sol_csv, file=file)


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
    """
    .DC analisia burutzen duen funtzioa. M, N eta Us matrizeak sortzen ditu
    uneko t balioarekin, sorgailu batean aplikatuz.

    Args
    ----------
    b : integer
        Gure zirkuituko adar kopurua.
    cir_el_extended : cir_el array luzatua.
    cir_val_extended : cir_val array luzatua.
    cir_ctr_extended : cir_ctr array luzatua.
    t : float
        Uneko balioa .DC azterketarako (tentsioa edo korrontea).
    sorgailua : str
        Baldintza aldakorra duen sorgailuaren izena.

    Returns
    -------
    lista : list
        M, N eta Us matrizeak bueltatzen ditu.
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
    j = np.where(np.char.lower(cir_el_extended) == np.char.lower(sorgailua))[0]
    Us[j] = t
    lista = []
    lista.append(M)
    lista.append(N)
    lista.append(Us)
    return lista


def op_tr(b, cir_el_extended, cir_val_extended, cir_ctr_extended, t):
    """
    .TRAN (transient) analisia burutzen duen funtzioa. M, N eta Us matrizeak
    sortzen ditu uneko t balioaren arabera, denboraren araberako portaera
    aztertzeko.

    Args
    ----------
    b : integer
        Gure zirkuituko adar kopurua.
    cir_el_extended : cir_el array luzatua.
    cir_val_extended : cir_val array luzatua.
    cir_ctr_extended : cir_ctr array luzatua.
    t : float
        Denbora balioa .TRAN azterketarako.

    Returns
    -------
    lista : list
        M, N eta Us matrizeak bueltatzen ditu.
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
            anp = float(cir_val_extended[i][0])
            maiz = float(cir_val_extended[i][1])
            fas = float(cir_val_extended[i][2])
            m = 2 * np.pi * maiz
            f = (np.pi / 180) * fas
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
    """

    Funtzio honek .cir fitxategi bateko balioak kargatzen ditu eta
    parametroak hartzen ditu.

    Args:
        filename: String-a fitxategiaren izenarekin

    Returns:
        A tuple with (sorgailua, hasiera, amaiera, pausua, a, pr, op).
    """
    sorgailua = None
    hasiera, amaiera, pausua = None, None, None
    a = None  # Analisi mota (.tr o .dc)
    pr = None  # Printa egin behar den adierazi
    op = None  # Operazio puntua egin behar den adierazi
    with open(filename, 'r') as file:
        lines = file.readlines()

        for line in lines:
            tokens = line.strip().split()

            if len(tokens) < 1:
                continue  # Linea hutsak kontuan ez
            
            if tokens[0].lower() == ".tr" or tokens[0].lower() == ".dc":
                a = tokens[0].lower()
                hasiera = float(tokens[5])  # Hasiera
                amaiera = float(tokens[6])  # Amaiera
                pausua = float(tokens[7])   # Denbora pausua
                sorgailua = str(tokens[8])

            if tokens[0].lower() == ".pr":
                pr = ".pr"

            if tokens[0].lower() == ".op":
                op = ".op"

    return sorgailua, hasiera, amaiera, pausua, a, pr, op
