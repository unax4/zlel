#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

.. module:: zlel_p3.py
    :synopsis: Zirkuitu erresistibo ez-linealak (diodoa eta transistorea)
               ebazteko gaitasuna gehituko diogu. Newton-Raphsonen metodoan eta
               baliokide diskretuan oinarrituko gara.

.. moduleauthor:: Aitor Lavin (aitorlavin02@gmail.com)
                  Unax Arregi (arregiunax@gmail.com)

"""

import numpy as np
import math

if __name__ == "zlel.zlel_p3":
    import zlel.zlel_p1 as zl1
    import zlel.zlel_p2 as zl2
else:
    import zlel_p1 as zl1
    import zlel_p2 as zl2


def diode_NR(I0, nD, Vdj):
    """ https://documentation.help/Sphinx/math.html
        Calculates the g and the I of a diode for a NR discrete equivalent
        Given,

        :math:`Id = I_0(e^{(\\frac{V_d}{nV_T})}-1)`

        The NR discrete equivalent will be,

        :math:`i_{j+1} + g v_{j+1} = I`

        where,

        :math:`g = -\\frac{I_0}{nV_T}e^{(\\frac{V_d}{nV_T})}`

        and

        :math:`I = I_0(e^{(\\frac{V_{dj}}{nV_T})}-1) + gV_{dj}`

    Args:
        | 0: Value of I0.
        | nD: Value of nD.
        | Vd: Value of Vd.

    Return:
        | gd: Conductance of the NR discrete equivalent for the diode.
        | Id: Current independent source of the NR discrete equivalent.

    """

    Vt = 8.6173324e-5*300*nD
    gd = -I0 * math.e ** (Vdj / Vt) / Vt
    Id = I0 * (math.e ** (Vdj/Vt) - 1) + gd*Vdj

    return gd, Id


def alpha(beta, Ies, Ics):
    """

   αF eta αR lortzeko funtzioa

   Args
   ----------
   beta : float bat transistorearen beta balioarekin.
   Ies : Float bat Ies korrontearen balioarena.
   Ics : Float bat Ics korrontearen balioarena.

   Returns
   -------
   alphaf : float bat alphaf balioarena
   alphar : float bat alphar balioarena

   """

    alphaf = beta/(1+beta)
    alphar = (Ies/Ics)*alphaf
    return alphaf, alphar


def transistoreko_balioak(Ies, Ics, Vbe, Vbc, alphaf, alphar):
    """

    transistorearen g11, g12, g21, g22, Ie, Ic balioak lortzeko funtzioa

    Args
    ----------
    Ies : Float bat Ies korrontearen balioarena.
    Ics : Float bat Ics korrontearen balioarena.
    Vbe : float balio bat Vbe tentsioaren balioarekin
    Vbc : float balio bat Vbc tentsioaren balioarekin
    alphaf : float bat alphaf balioarena
    alphar : float bat alphar balioarena

    Returns
    -------
    g11 : float bat g11 balioarekin
    g12 : float bat g11 balioarekin
    g21 : float bat g21 balioarekin
    g22 : float bat g22 balioarekin
    Ie : float bat transistorearen igorleko korrontearen balioarekin
    Ic : float bat transistorearen kolektoreko korrontearen balioarekin

    """

    Vt = 8.6173324e-5*300
    g11 = -Ies/Vt*math.e**(Vbe/Vt)
    g22 = -Ics/Vt*math.e**(Vbc/Vt)
    g12 = -alphar*g22
    g21 = -alphaf*g11
    Ie = g11*Vbe + g12*Vbc + Ies*(math.e**(Vbe/Vt)-1) -\
        alphar*Ics*(math.e**(Vbc/Vt)-1)
    Ic = g21*Vbe + g22*Vbc - alphaf*Ies*(math.e**(Vbe/Vt)-1) +\
        Ics*(math.e**(Vbc/Vt)-1)

    return g11, g12, g21, g22, Ie, Ic


def elementu_ezlinealak(b, cir_el_extended):
    """

    Aztertu ea elementu ez-linealak dauden ala ez, eta zein adarretan dauden.

    Args
    ----------
    b : integer-a adar kopuruarena.
    cir_el_extended : cir_el array luzatua.

    Returns
    -------
    list
        zerrenda bat bueltatuko du beste bi zerrendekin bere barruan. 1
        diodoekin eta bestea transistoreekin.

    """

    diodo_lista = []
    transistore_lista = []
    i = 0
    while i < b:
        if cir_el_extended[i][0].lower() == "d":
            diodo_lista.append(i)
            i += 1
        elif cir_el_extended[i][0].lower() == "q":
            transistore_lista.append(i)
            i += 2
        else:
            i += 1
    return [diodo_lista, transistore_lista]


def Newton_Raphson(M, N, Us, A, At, n, b, cir_val_extended,
                   diodo_lista, transistore_lista):
    """
    Elementu ez-linealak baldin badaude, Newton-Raphsonen metodoa aplikatuko
    du funtzio honek.

    Args
    ----------
    M : np array-a M matrizearena.
    N : np array-a N matrizearena.
    Us : np array-a Us matrizearena.
    A : np array-a M matrizearen matrize murriztuarena
    At : np_array-a A matrizearen irauliarena.
    n : Integer-a nodo kopuruarekin.
    b : Integer-a adar kopuruarekin
    cir_val_extended : cir_val array luzatua
    diodo_lista : Zerrenda bat zirkuituko diodoena
    transistore_lista : Zerrenda bat zirkuituko transistoreena.

    Returns
    -------
    sol : Newton-Rhapson metodoaren emaitza.

    """

    i = 0
    Vd_lista = []
    k_lista = []
    for h in range(0, b):
        k_lista.append(0)
        if h in diodo_lista:
            Vd_lista.append(0.7)
            k_lista.append(0)
        elif h in transistore_lista:
            Vd_lista.append(0.6)
            Vd_lista.append(0.6)
            k_lista.append(0)
            k_lista.append(0)
            h += 1
        else:
            Vd_lista.append(0)
            k_lista.append(0)
    # Newton Ramphson metodoa elementu ezberdinekin
    while i < 100:
        i += 1
        for j in diodo_lista:
            gd, Id = diode_NR(float(cir_val_extended[j][0]),
                              float(cir_val_extended[j][1]), Vd_lista[j])
            M[j][j] = gd
            N[j][j] = 1
            Us[j] = Id
        for j in transistore_lista:
            alphaf, alphar = alpha(float(cir_val_extended[j][2]),
                                   float(cir_val_extended[j][0]),
                                   float(cir_val_extended[j][1]))
            g11, g12, g21, g22, Ie, Ic = transistoreko_balioak(
                float(cir_val_extended[j][0]), float(cir_val_extended[j][1]),
                float(Vd_lista[j]), Vd_lista[j + 1], alphaf, alphar)
            M[j][j], M[j][j+1], M[j+1][j], M[j+1][j+1] = g11, g12, g21, g22
            N[j][j], N[j+1][j+1] = 1, 1
            Us[j], Us[j+1] = Ie, Ic
        T = zl2.get_T(A, At, M, N, Us, b, n)
        U = zl2.get_U(A, b, n, Us)
        sol = np.linalg.solve(T, U)
        for j in diodo_lista:
            if k_lista[j] == 0:
                Vd = sol[n-1+j]
                if abs(Vd - Vd_lista[j]) < 10 ** (-5):
                    k_lista[j] = 1
                Vd_lista[j] = float(Vd)
        for j in transistore_lista:
            if k_lista[j] == 0:
                Vd1, Vd2 = sol[n-1+j], sol[n+j]
                if abs(Vd1 - Vd_lista[j]) < 10 ** (-5) and\
                        abs(Vd2 - Vd_lista[j+1]) < 10 ** (-5):
                    k_lista[j] = 1
                Vd_lista[j], Vd_lista[j+1] = float(Vd1), float(Vd2)
        if all(k_lista[j] == 1 for j in diodo_lista + transistore_lista):
            return sol


def save_as_csv_p3(b, n, filename, cir_el_extended, cir_val_extended,
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
    header = zl2.build_csv_header("t", b, n)
    with open(filename, 'w') as file:
        print(header, file=file)
        t = hasiera
        num_steps = int(round((amaiera - hasiera) / pausua)) + 1
        for i in range(num_steps):
            t = hasiera + i * pausua
            if a == ".tr":
                m = zl2.op_tr(b, cir_el_extended, cir_val_extended,
                              cir_ctr_extended, t)
            elif a == ".dc":
                m = zl2.op_dc(b, cir_el_extended, cir_val_extended,
                              cir_ctr_extended, t, sorgailua)
            M = m[0]
            N = m[1]
            Us = m[2]
            zl1.erroreak(cir_nd_extended, cir_el_extended, cir_ctr_extended,
                         cir_val_extended, incidence_matrix, nodes, b)
            diodo_lista, transistore_lista = elementu_ezlinealak(
                b, cir_el_extended)
            sol = Newton_Raphson(M, N, Us, A, At, n, b, cir_val_extended,
                                 diodo_lista, transistore_lista)
            sol = np.insert(sol, 0, t)            
            sol_csv = ','.join(['%.9f' % num for num in sol])
            print(sol_csv, file=file)
