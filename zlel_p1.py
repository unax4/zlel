#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: zlel_main.py
    :synopsis:

.. moduleauthor:: YOUR NAME AND E-MAIL


"""

import numpy as np
import sys


def cir_parser(filename):
    """
        This function takes a .cir test circuit and parse it into
        4 matices.
        If the file has not the proper dimensions it warns and exit.

    Args:
        filename: string with the name of the file

    Returns:
        | cir_el: np array of strings with the elements to parse. size(1,b)
        | cir_nd: np array with the nodes to the circuit. size(b,4)
        | cir_val: np array with the values of the elements. size(b,3)
        | cir_ctrl: np array of strings with the element which branch
        | controls the controlled sources. size(1,b)

    Rises:
        SystemExit

    """
    try:
        cir = np.array(np.loadtxt(filename, dtype=str))
    except ValueError:
        sys.exit("File corrupted: .cir size is incorrect.")

    # numpy usefull exmaples
    # print("================ cir ==========")
    # print(cir)
    # print("\n======== a = np.array (cir[:,1], dtype = int) ==========")
    # a = np.array(cir[:, 1], dtype=int)
    # print(a)
    # print("\n======== a = np.append(a,300) ==========")
    # a = np.append(a, 300)
    # print(a)
    # print("\n======== b = a[a > 3] ==========")
    # b = a[a > 3]
    # print(b)
    # print("\n======== c = np.unique(a) ==========")
    # c = np.unique(a)
    # print(c)
    # print("\n======== d = np.flatnonzero(a != 0) ==========")
    # d = np.flatnonzero(a != 0)
    # print(d)
    # print("\n======== e = np.flatnonzero(a == 0) ==========")
    # e = np.flatnonzero(a == 0)
    # print(e)
    # print("\n======== f = np.array(cir[:, 1:2]) ==========")
    # f = np.array(cir[:, 1:2])
    # print(f)
    # print("\n======== g = np.array(cir[2:4, :]) ==========")
    # g = np.array(cir[2:4, :])
    # print(g)
    # print("\n======== h = np.empty([0], dtype=int) ==========")
    # h = np.empty([0], dtype=int)
    # print(h)
    # print("\n======== i = np.append(h, 1) ==========")
    # i = np.append(h, 1)
    # print(i)
    # print("\n======== i[0] = 2 ==========")
    # i[0] = 2
    # print(i)
    # print("\n======== j = np.empty([0], dtype=str ==========")
    # j = np.empty([0], dtype=str)
    # print(j)
    # print("\n======== k = np.append(j, \"123456\") ==========")
    # k = np.append(j, "123456")
    # print(k)
    # print("\n======== k[0] = \"987654321\" ==========")
    # k[0] = "987654321"
    # print(k)
    ''' https://www.geeksforgeeks.org/modify-numpy-array-to-store-an-arbitrary-
    length-string/
    The dtype of any numpy array containing string values is the maximum
    length of any string present in the array. Once set, it will only be able
    to store new string having length not more than the maximum length at the
    time of the creation. If we try to reassign some another string value
    having length greater than the maximum length of the existing elements,
    it simply discards all the values beyond the maximum length.'''

    cir_el = np.array(cir[:-1, 0], dtype=str)
    cir_nd = np.array(cir[:-1, 1:5], dtype=int)
    cir_val = np.array(cir[:-1, 5:8])
    cir_ctr = np.array(cir[:-1, 8], dtype=str)
    cir_el_ona = []
    cir_nd_ona = []
    cir_val_ona = []
    cir_ctr_ona = []
    for el, nd, val, ctr in zip(cir_el, cir_nd, cir_val, cir_ctr):
        if el[0].lower() == ".":
            continue
        else:
            cir_el_ona.append(el)
            cir_nd_ona.append(nd)
            cir_val_ona.append(val)
            cir_ctr_ona.append(ctr)
    return (np.array(cir_el_ona), np.array(cir_nd_ona),
            np.array(cir_val_ona), np.array(cir_ctr_ona))


def modify_function(cir_el, cir_nd, cir_val, cir_ctr):
    """
    Funtzio honek aurretik sortutako cir_el, cir_nd... hartzen
    ditu eta array-a luzatzen du q edo a-ren arabera.

    Parameters
    ----------
    cir_el : np array of strings with the elements to parse. size(1,b)
    cir_nd : np array with the nodes to the circuit. size(b,4)
    cir_val : np array with the values of the elements. size(b,3)
    cir_ctr : np array of strings with the element which branch
    controls the controlled sources. size(1,b)

    Returns
    -------
    cir_el_ext : cir_el array luzatua
    cir_nd_ext : cir_nd array luzatua
    cir_val_ext : cir_val array luzatua
    cir_ctr_ext : cir_ctr array luzatua

    """

    cir_el_ext = []
    cir_nd_ext = []
    cir_val_ext = []
    cir_ctr_ext = []

    for el, nd, val, ctrl in zip(cir_el, cir_nd, cir_val, cir_ctr):
        if el[0].lower() == "q":
            cir_el_ext.extend([el + "_be", el + "_bc"])
            cir_nd_ext.extend([[nd[1], nd[2], 0, 0], [nd[1], nd[0], 0, 0]])
            cir_val_ext.extend([val, val])
            cir_ctr_ext.extend([ctrl, ctrl])
        elif el[0].lower() == "a":
            cir_el_ext.extend([el + "_in", el + "_ou"])
            cir_nd_ext.extend([[nd[0], nd[1], 0, 0], [nd[2], nd[3], 0, 0]])
            cir_val_ext.extend([val, val])
            cir_ctr_ext.extend([ctrl, ctrl])
        elif el[0].lower() == ".":
            continue
        else:
            cir_el_ext.append(el)
            cir_nd_ext.append(nd)
            cir_val_ext.append(val)
            cir_ctr_ext.append(ctrl)

    return cir_el_ext, cir_nd_ext, cir_val_ext, cir_ctr_ext


def get_branches(cir_el):
    """
    Gure zirkuituaren adar kopurua ematen du.

    Parameters
    ----------
    cir_el : np array-a gure zirkuituko elementuen lista

    Returns
    -------
    b : integer
        Gure zirkuituko adar kopurua.

    """
    b = 0
    for x in cir_el:
        if x[0].lower() == "q" or x[0].lower() == "a":
            b += 2
        else:
            b += 1
    return b


def get_nodes(cir_nd):
    """
    Zirkuituko nodoen lista bueltatzen du.

    Parameters
    ----------
    cir_nd : np array with the nodes to the circuit. size(b,4)

    Returns
    -------
    nodes : np array-a zirkuituko nodoekin.

    """
    all_nodes = np.unique(cir_nd)
    valid_nodes = set(all_nodes)
    if 0 in cir_nd[:, :2]:
        pass
    else:
        valid_nodes.discard(0)
    return np.array(sorted(valid_nodes))


def get_number_of_nodes(nodes):
    """
    Nodo kopurua ematen du.

    Parameters
    ----------
    nodes : np array-a zirkuituko nodoekin.

    Returns
    -------
    n : Integer-a nodo kopuruarekin

    """
    n = nodes.size
    return (n)


def get_elements(cir_el):
    """
    Zirkuituko elementu kopurua ematen digu.

    Parameters
    ----------
    cir_el : np array-a gure zirkuituko elementuen lista

    Returns
    -------
    el_num : Integer-a zirkuituko elementu kopuruarena.

    """
    el_num = cir_el.size
    return (el_num)


def get_incidence_matrix(n, b, nodes, cir_nd_extended):
    """
    Funtzio honek zirkuituaren intzidentzia matrizea sortzen du.

    Parameters
    ----------
    n : Integer-a nodo kopuruarekin
    b : Integer-a zirkuituko elementu kopuruarena.
    nodes : np array-a zirkuituko nodoekin.
    cir_nd_extended : cir_nd array luzatua

    Returns
    -------
    incidence_matrix : np array-a zirkuituaren intzidentzia matrizearena.

    """
    incidence_matrix = np.zeros((n, b), dtype=int)
    for i in range(b):
        node_pos = np.where(nodes == cir_nd_extended[i, 0])[0][0]
        incidence_matrix[node_pos, i] = 1
        node_pos = np.where(nodes == cir_nd_extended[i, 1])[0][0]
        incidence_matrix[node_pos, i] = -1
    return incidence_matrix


"INTZIDENTZIA MATRIZE MURRIZTUA"


def get_incidence_murriztua(incidence_matrix, n):
    """
   Intzidentzia matrize murriztua lortzen du.

   Parameters
   ----------
   incidence_matrix : np array-a zirkuituaren intzidentzia matrizearena.
   n : Integer-a nodo kopuruarekin

   Returns
   -------
   np array-a intzidentzia matrize murriztuarekin.

   """
    return np.array(incidence_matrix[1:n, :])


def branchZerrenda(nodoa, cir_el_extended, cir_nd_extended):
    """
    Zirkuituaren adarren zerrenda ematen du.

    Parameters
    ----------
    nodoa : Integer-a nodoa adierazten duena
    cir_el_extended : cir_el array luzatua.
    cir_nd_extended : cir_nd array luzatua

    Returns
    -------
    lista : Zerrenda bat bueltatzen du zirkuituko adarrekin.

    """
    lista = []
    for i in range(0, np.size(cir_el_extended)):
        for x in cir_nd_extended[i]:
            if x == nodoa:
                lista.append(cir_el_extended[i])
    return lista


def erroreak(cir_nd_extended, cir_el_extended, cir_ctr_extended,
             cir_val_extended, incidence_matrix, nodes, b):
    """
    5 errore posible ezberdinen artetik zein den bueltatuko duen funtzioa da.

    Parameters
    ----------
    cir_nd_extended : cir_nd array luzatua
    cir_el_extended : cir_el array luzatua
    cir_ctr_extended : cir_ctr array luzatua
    cir_val_extended : cir_val array luzatua
    incidence_matrix : np array-a zirkuituaren intzidentzia matrizearena.
    nodes : np array-a zirkuituko nodoekin.
    b : Integer-a zirkuituko elementu kopuruarena.

    Returns
    -------
    String bat errorea azalduz.

    """
    # Ez dago erreferentzia nodorik
    i = 0
    for x in nodes:
        if x == 0:
            i = 1
    if i == 0:
        sys.exit("Reference node \"0\" is not defined in the circuit.")

    # Tentsio iturriak paraleloan konektatuta.
    cir_val_extended = np.array(cir_val_extended, dtype=np.float64)
    v = 0
    hiztegia = {}
    for x in cir_el_extended:
        if x[0].lower() in ("v", "e", "b"):
            hiztegia[v] = cir_nd_extended[v]
        v += 1
    for x in hiztegia.keys():
        for y in hiztegia.keys():
            if x != y:
                if hiztegia[x][0] == hiztegia[y][0] and \
                        hiztegia[x][1] == hiztegia[y][1]:
                    if cir_val_extended[x][0] != cir_val_extended[y][0]:
                        sys.exit("Parallel V sources at branches " + str(x) +
                                 " and " + str(y) + ".")
                if hiztegia[x][0] == hiztegia[y][1] and \
                        hiztegia[x][0] == hiztegia[y][1]:
                    if cir_val_extended[x][0] != -cir_val_extended[y][0]:
                        sys.exit("Parallel V sources at branches " + str(x) +
                                 " and " + str(y) + ".")

    # Korronte iturriak seriean konektatuta.
    v = 0
    hiztegia = {}
    nodoZerrenda = []
    for x in cir_el_extended:
        if x[0].lower() in ("i", "g", "y"):
            hiztegia[v] = cir_nd_extended[v]
        v += 1
    for x in hiztegia.keys():
        for y in hiztegia.keys():
            if x != y:
                if hiztegia[x][0] == hiztegia[y][0] and \
                        hiztegia[x][1] != hiztegia[y][1]:
                    if hiztegia[x][0] not in nodoZerrenda:
                        nodoZerrenda.append(hiztegia[x][0])
                if hiztegia[x][0] == hiztegia[y][1] and \
                        hiztegia[x][1] != hiztegia[y][0]:
                    if hiztegia[x][0] not in nodoZerrenda:
                        nodoZerrenda.append(hiztegia[x][0])
                if hiztegia[x][1] == hiztegia[y][0] and\
                        hiztegia[x][0] != hiztegia[y][1]:
                    if hiztegia[x][0] not in nodoZerrenda:
                        nodoZerrenda.append(hiztegia[x][1])
                if hiztegia[x][1] == hiztegia[y][1] and \
                        hiztegia[x][0] != hiztegia[y][0]:
                    if hiztegia[x][0] not in nodoZerrenda:
                        nodoZerrenda.append(hiztegia[x][1])
    for x in nodoZerrenda:
        bZerrenda = branchZerrenda(x, cir_el_extended, cir_nd_extended)
        k = 0
        s = 0
        for y in bZerrenda:
            if y[0].lower() not in ("i", "g", "y"):
                k = 1
                break
        if k == 0:
            s = 0
            for i in range(0, b):
                s += cir_val_extended[i][0]*incidence_matrix[x][i]
            if s != 0:
                sys.exit("I sources in series at node " + str(x) + ".")

    # Nodoa konektatu gabe.
    j = 0
    for x in incidence_matrix:
        if np.size(np.flatnonzero(x != 0)) < 2:
            h = nodes[j]
            sys.exit("Node " + str(h) + " is floating.")
        j += 1


def print_cir_info(cir_el, cir_nd, b, n, nodes, el_num):
    """ Prints the info of the circuit:
        |     1.- Elements info
        |     2.- Node info
        |     3.- Branch info
        |     4.- Variable info
    Args:
        | cir_el: reshaped cir_el
        | cir_nd: reshaped cir_nd. Now it will be a(b,2) matrix
        | b: # of branches
        | n: # number of nodes
        | nodes: an array with the circuit nodes sorted
        | el_num:  the # of elements

    """
    # Element info
    print(str(el_num) + ' Elements')
    # Node info
    print(str(n) + ' Different nodes: ' +
          str(nodes))
    # Branch info
    print("\n" + str(b) + " Branches: ")

    for i in range(1, b+1):
        indent = 12  # Number of blanks for indent
        string = ("\t" + str(i) + ". branch:\t" +
                  str(cir_el[i-1]) + "i".rjust(indent - len(cir_el[i-1])) +
                  str(i) + "v".rjust(indent - len(str(i))) + str(i) +
                  " = e" + str(cir_nd[i-1, 0]) +
                  " - e" + str(cir_nd[i-1, 1]))
        print(string)

    # Variable info
    print("\n" + str(2*b + (n-1)) + " variables: ")
    # print all the nodes but the first(0 because is sorted)
    for i in nodes[1:]:
        print("e"+str(i)+", ", end="", flush=True)
    for i in range(b):
        print("i"+str(i+1)+", ", end="", flush=True)
    # print all the branches but the last to close it properly
    # It works because the minuimum amount of branches in a circuit must be 2.
    for i in range(b-1):
        print("v"+str(i+1)+", ", end="", flush=True)
    print("v"+str(b))

    # IT IS RECOMMENDED TO USE THIS FUNCTION WITH NO MODIFICATION.


"""
https://stackoverflow.com/questions/419163/what-does-if-name-main-do
https://stackoverflow.com/questions/19747371/
python-exit-commands-why-so-many-and-when-should-each-be-used
"""
if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "C:/Users/usuario/Desktop/3maila/Zlel/Aitor_Lavin.Unax_Arregi/cirs/all/0_zlel_OPAMP.cir"

    cir_el, cir_nd, cir_val, cir_ctr = cir_parser(filename)
    modified_results = modify_function(cir_el, cir_nd, cir_val, cir_ctr)
    cir_el_extended = np.array(modified_results[0])
    cir_nd_extended = np.array(modified_results[1])
    cir_ctr_extended = np.array(modified_results[3])
    cir_val_extended = np.array(modified_results[2])
    nodes = get_nodes(cir_nd)
    b = get_branches(cir_el)
    n = get_number_of_nodes(nodes)
    el_num = get_elements(cir_el)
    incidence_matrix = get_incidence_matrix(n, b, nodes, cir_nd_extended)
    erroreak(cir_nd_extended, cir_el_extended, cir_ctr_extended,
             cir_val_extended, incidence_matrix, nodes, b)
    print_cir_info(cir_el_extended, cir_nd_extended, b, n, nodes, el_num)
    print("\nIncidence Matrix:")
    print(incidence_matrix)
    erroreak(cir_nd_extended, cir_el_extended, cir_ctr_extended,
             cir_val_extended, incidence_matrix, nodes, b)
