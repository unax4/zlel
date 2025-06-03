#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: zlel_p4.py
    :synopsis: Euler Backward erabiliz, zirkuitu dinamikoen analisia egiten
               duen modulua.

.. moduleauthor:: Aitor Lavin (aitorlavin02@gmail.com)
                  Unax Arregi (arregiunax@gmail.com)
"""

import numpy as np

if __name__ == "zlel.zlel_p4":
    import zlel.zlel_p1 as zl1
    import zlel.zlel_p2 as zl2
    import zlel.zlel_p3 as zl3
else:
    import zlel_p1 as zl1
    import zlel_p2 as zl2
    import zlel_p3 as zl3


def capacitor_NR(C, Vc_prev, h):
    """

    Kapazitore baten NR baliokidearen parametroak kalkulatzen ditu,
    Euler Backward erabiliz.

    Args:
        C: Kondentsadorearen kapazitantzia
        Vc_prev: Kondentsadorearen tentsioa aurreko etapan
        h: denbora

    Returns:
        gc: Konduktzantzia baliokidea
        Ic: Korronte iturri baliokidea
    """
    gc = C / h
    Ic = gc * Vc_prev
    return gc, Ic


def inductor_NR(L, Il_prev, h):
    """

    Haril baten NR baliokidearen parametroak kalkulatzen ditu,
    Euler Backward erabiliz.

    Args:
        L: Harilaren induktantzia
        Il_prev: Harilaren korrontea aurreko etapan
        h: denbora

    Returns:
        gl: Konduktzantzia baliokidea
        Il: Korronte iturri baliokidea
    """
    gl = h / L
    Vl = Il_prev
    return gl, Vl


def elementu_dinamikoak(b, cir_el_extended):
    """

    Zirkuituko elemntu dinamikoa identifikatzen ditu (kondentsadore
    eta induktoreak).

    Args:
        b: Adar kopurua
        cir_el_extended: Cir_el array luzatua

    Returns:
        Bi subzerrenda dituen zerrenda: [kondentsadoreak, induktoreak]
    """
    capacitor_lista = []
    inductor_lista = []

    for i in range(b):
        if cir_el_extended[i][0].lower() == "c":
            capacitor_lista.append(i)
        elif cir_el_extended[i][0].lower() == "l":
            inductor_lista.append(i)

    return [capacitor_lista, inductor_lista]


def dynamic_NR(b, cir_el_extended, cir_val_extended, h, Vc_prev, Il_prev,
               capacitor_lista, inductor_lista, M, N, Us):
    """

    M, N, Us matrizeak eguneratzen ditu elementu dinamikoentzat,
    Euler Backward erabiliz.

    Args:
        b: Adar kopurua
        cir_el_extended: Cir_el array luzatua
        cir_val_extended: Cir_val array luzatua
        h: denbora
        Vc_prev: Kondentsadorearen tentsioa aurreko etapan
        Il_prev: Harilaren korrontea aurreko etapan
        capacitor_lista: Kondentsadoreen zerrenda
        inductor_lista: Harilen zerrenda
        M: np array-a M matrizearena.
        N: np array-a N matrizearena.
        Us: np array-a Us matrizearena.

    Returns:
        M: np array-a M matrizearena.
        N: np array-a N matrizearena.
        Us: np array-a Us matrizearena.
        Vc_actual: Kondentsadorearen tentsioa etapa honetan
        Il_actual: Harilaren korrontea etapa honetan
    """
    Vc_actual = np.zeros(b)
    Il_actual = np.zeros(b)

    for j in capacitor_lista:
        C = float(cir_val_extended[j][0])
        if len(cir_val_extended[j]) > 1:
            Vc0 = float(cir_val_extended[j][1])
        else:
            Vc0 = 0.0

        if np.all(Vc_prev == 0) and Vc0 != 0:
            Vc_prev[j] = Vc0

        gc, Ic = capacitor_NR(C, Vc_prev[j], h)
        M[j, j] = gc
        N[j, j] = -1
        Us[j, 0] = Ic
        Vc_actual[j] = Vc_prev[j]

    for j in inductor_lista:
        L = float(cir_val_extended[j][0])
        if len(cir_val_extended[j]) > 1:
            Il0 = float(cir_val_extended[j][1])
        else:
            Il0 = 0.0
        if Il_prev[j] == 0 and Il0 != 0:
            Il_prev[j] = Il0

        gl, Vl = inductor_NR(L, Il_prev[j], h)
        M[j, j] = -gl
        N[j, j] = 1
        Us[j, 0] = Vl
        Il_actual[j] = Il_prev[j]

    return M, N, Us, Vc_actual, Il_actual


def save_as_csv_p4(b, n, filename, cir_el_extended, cir_val_extended,
                   cir_ctr_extended, A, At, cir_nd, cir_nd_extended,
                   incidence_matrix, nodes, a, hasiera, amaiera, pausua,
                   sorgailua):
    """
    This function generates a csv file with the name filename.
        First it will save a header and then, it loops and save a line in
        csv format into the file.

    Args:
        b: Adar kopurua
        n: Nodo kopurua
        filename: Irteerako fitxategi izena
        A: Inzidentzia matrize murriztua
        At: A-ren iraulia
        incidence_matrix: Inzidentzia matrizea
        nodes: Nodo lista
        a: Analisi mota (.tr o .dc)
        hasiera: Hasiera denbora
        amaiera: Amaiera denbora
        pausua: Denbora pausua
        sorgailua: Kontrol iturria
    """
    header = zl2.build_csv_header("t", b, n)
    capacitor_lista, inductor_lista = elementu_dinamikoak(b, cir_el_extended)
    diodo_lista, transistore_lista = zl3.elementu_ezlinealak(b,
                                                             cir_el_extended)

    Vc_prev = np.zeros(b)
    Il_prev = np.zeros(b)

    num_steps = int(np.ceil((amaiera - hasiera) / pausua)) + 1
    with open(filename, 'w') as file:
        print(header, file=file)

        for step in range(num_steps):
            t = hasiera + step * pausua
            if abs(t - amaiera) < 1e-10 or t <= amaiera:
                if a == ".tr":
                    m = zl2.op_tr(b, cir_el_extended, cir_val_extended,
                                  cir_ctr_extended, t)
                elif a == ".dc":
                    m = zl2.op_dc(b, cir_el_extended, cir_val_extended,
                                  cir_ctr_extended, t, sorgailua)

                M, N, Us = m[0], m[1], m[2]

                zl1.erroreak(cir_nd_extended, cir_el_extended,
                             cir_ctr_extended, cir_val_extended,
                             incidence_matrix, nodes, b)

                M, N, Us, Vc_actual, Il_actual = dynamic_NR(
                    b, cir_el_extended, cir_val_extended, pausua,
                    Vc_prev, Il_prev, capacitor_lista,
                    inductor_lista, M, N, Us)

                sol = zl3.Newton_Raphson(M, N, Us, A, At,
                                         n, b, cir_val_extended,
                                         diodo_lista, transistore_lista)

                for j in capacitor_lista:
                    Vc_prev[j] = float(sol[n-1 + j])

                for j in inductor_lista:
                    Il_prev[j] = float(sol[n-1 + b + j])

                sol = np.insert(sol, 0, t)
                sol_csv = ','.join(['%.9f' % num for num in sol])
                print(sol_csv, file=file)
