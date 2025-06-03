#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: zlel_p5.py
    :synopsis: Modulua denborarekiko aldatzen den erresistentzia (RV) elementuak
               kudeatzeko, .tr simulazioetan. Erresistentzia R0 balioa du t <= Time
               eta RF balioa t > Time.

.. moduleauthor:: Aitor Lavin (aitorlavin02@gmail.com)
                  Unax Arregi (arregiunax@gmail.com)
"""

import numpy as np

if __name__ == "zlel.zlel_p5":
    import zlel.zlel_p1 as zl1
    import zlel.zlel_p2 as zl2
    import zlel.zlel_p3 as zl3
    import zlel.zlel_p4 as zl4
else:
    import zlel_p1 as zl1
    import zlel_p2 as zl2
    import zlel_p3 as zl3
    import zlel_p4 as zl4


def time_varying_resistor(b, cir_el_extended):
    """
    Identifikatu zirkuituko denborarekiko aldatzen diren erresistentziak (RV).

    Args:
        b: Adar kopurua
        cir_el_extended: Cir_el array luzatua

    Returns:
        rv_lista: Denborarekiko aldatzen diren erresistentzien indizeen zerrenda
    """
    rv_lista = []
    for i in range(b):
        if cir_el_extended[i][0].lower() == "rv":
            rv_lista.append(i)
    return rv_lista


def op_tr_p5(b, cir_el_extended, cir_val_extended, cir_ctr_extended, t):
    """
    .TRAN (transient) analisia burutzen duen funtzioa, RV elementuak kontuan hartuz.
    M, N eta Us matrizeak sortzen ditu uneko t balioaren arabera, RV elementuen
    erresistentzia R0 edo RF balioetara egokituz t > Time den ala ez.

    Args:
        b: Gure zirkuituko adar kopurua.
        cir_el_extended: Cir_el array luzatua.
        cir_val_extended: Cir_val array luzatua.
        cir_ctr_extended: Cir_ctr array luzatua.
        t: Denbora balioa .TRAN azterketarako.

    Returns:
        lista: M, N eta Us matrizeak bueltatzen ditu.
    """
    M = np.zeros((b, b), dtype=float)
    N = np.zeros((b, b), dtype=float)
    Us = np.zeros((b, 1), dtype=float)
    rv_lista = time_varying_resistor(b, cir_el_extended)

    for i in range(b):
        if cir_el_extended[i][0].lower() == "r":
            M[i][i] = 1
            N[i][i] = -float(cir_val_extended[i][0])
        elif cir_el_extended[i][0].lower() == "rv":
            R0 = float(cir_val_extended[i][0])
            RF = float(cir_val_extended[i][1])
            Time = float(cir_val_extended[i][2])
            M[i][i] = 1
            N[i][i] = -RF if t > Time else -R0
        elif cir_el_extended[i][0].lower() == "v":
            M[i][i] = 1
            Us[i] = float(cir_val_extended[i][0])
        elif cir_el_extended[i][0].lower() == "i":
            N[i][i] = 1
            Us[i] = float(cir_val_extended[i][0])
        elif cir_el_extended[i][0].lower() == "a":
            if "in" in cir_el_extended[i]:
                M[i][i] = 1
                N[i+1][i] = 1
        elif cir_el_extended[i][0].lower() == "e":
            j = np.where(np.char.lower(cir_el_extended) == np.char.lower(
                cir_ctr_extended[i]))[0]
            M[i][i] = 1
            M[i][j] = -float(cir_val_extended[i][0])
        elif cir_el_extended[i][0].lower() == "g":
            j = np.where(np.char.lower(cir_el_extended) == np.char.lower(
                cir_ctr_extended[i]))[0]
            N[i][i] = 1
            M[i][j] = -float(cir_val_extended[i][0])
        elif cir_el_extended[i][0].lower() == "f":
            j = zl2.getPosition(cir_ctr_extended[i], cir_el_extended)
            N[i][i] = 1
            N[i][j] = -float(cir_val_extended[i][0])
        elif cir_el_extended[i][0].lower() == "h":
            j = zl2.getPosition(cir_ctr_extended[i], cir_el_extended)
            M[i][i] = 1
            N[i][j] = -float(cir_val_extended[i][0])
        elif cir_el_extended[i][0].lower() == "b":
            anp = float(cir_val_extended[i][0])
            maiz = float(cir_val_extended[i][1])
            fas = float(cir_val_extended[i][2])
            m = 2 * np.pi * maiz
            f = (np.pi / 180) * fas
            M[i, i] = 1
            Us[i] = anp * np.sin(m * t + f)
        elif cir_el_extended[i][0].lower() == "y":
            anp = float(cir_val_extended[i][0])
            maiz = float(cir_val_extended[i][1])
            fas = float(cir_val_extended[i][2])
            m = 2 * np.pi * maiz
            f = float((np.pi/180) * fas)
            N[i, i] = 1
            Us[i] = anp * np.sin(m * t + f)

    return [M, N, Us]


def save_as_csv_p5(b, n, filename, cir_el_extended, cir_val_extended,
                   cir_ctr_extended, A, At, cir_nd, cir_nd_extended,
                   incidence_matrix, nodes, a, hasiera, amaiera, pausua,
                   sorgailua):
    """
    CSV fitxategia sortzen du .tr simulazioetarako, RV elementuak kontuan hartuz.
    Lehenik, goiburua gordetzen du eta, ondoren, lerroak CSV formatuan.

    Args:
        b: Adar kopurua
        n: Nodo kopurua
        filename: Irteerako fitxategi izena
        cir_el_extended: Cir_el array luzatua
        cir_val_extended: Cir_val array luzatua
        cir_ctr_extended: Cir_ctr array luzatua
        A: Inzidentzia matrize murriztua
        At: A-ren iraulia
        cir_nd: Nodoen array-a
        cir_nd_extended: Nodoen array luzatua
        incidence_matrix: Inzidentzia matrizea
        nodes: Nodo lista
        a: Analisi mota (.tr)
        hasiera: Hasiera denbora
        amaiera: Amaiera denbora
        pausua: Denbora pausua
        sorgailua: Kontrol iturria
    """
    header = zl2.build_csv_header("t", b, n)
    capacitor_lista, inductor_lista = zl4.elementu_dinamikoak(b, cir_el_extended)
    diodo_lista, transistore_lista = zl3.elementu_ezlinealak(b, cir_el_extended)
    rv_lista = time_varying_resistor(b, cir_el_extended)

    Vc_prev = np.zeros(b)
    Il_prev = np.zeros(b)

    num_steps = int(np.ceil((amaiera - hasiera) / pausua)) + 1
    with open(filename, 'w') as file:
        print(header, file=file)

        for step in range(num_steps):
            t = hasiera + step * pausua
            if abs(t - amaiera) < 1e-10 or t <= amaiera:
                m = op_tr_p5(b, cir_el_extended, cir_val_extended,
                             cir_ctr_extended, t)
                M, N, Us = m[0], m[1], m[2]

                zl1.erroreak(cir_nd_extended, cir_el_extended,
                             cir_ctr_extended, cir_val_extended,
                             incidence_matrix, nodes, b)

                if capacitor_lista or inductor_lista:
                    M, N, Us, Vc_actual, Il_actual = zl4.dynamic_NR(
                        b, cir_el_extended, cir_val_extended, pausua,
                        Vc_prev, Il_prev, capacitor_lista,
                        inductor_lista, M, N, Us)
                else:
                    Vc_actual, Il_actual = np.zeros(b), np.zeros(b)

                if diodo_lista or transistore_lista:
                    sol = zl3.Newton_Raphson(M, N, Us, A, At, n, b,
                                             cir_val_extended, diodo_lista,
                                             transistore_lista)
                else:
                    T = zl2.get_T(A, At, M, N, Us, b, n)
                    U = zl2.get_U(A, b, n, Us)
                    zl2.SingularMatrix(T, U)
                    sol = np.linalg.solve(T, U)

                for j in capacitor_lista:
                    Vc_prev[j] = float(sol[n-1 + j])
                for j in inductor_lista:
                    Il_prev[j] = float(sol[n-1 + b + j])

                sol = np.insert(sol, 0, t)
                sol_csv = ','.join(['%.9f' % num for num in sol])
                print(sol_csv, file=file)

if __name__ == "__main__":
    import sys
    import os
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "D://Universidad//Zlel//Praktikak//arregi_lavin//zlel//adb.cir"

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

    zl1.erroreak(cir_nd_extended, cir_el_extended, cir_ctr_extended,
                 cir_val_extended, incidence_matrix, nodes, b)

    A = zl1.get_incidence_murriztua(incidence_matrix, n)
    At = zl2.get_iraulia(A)

    diodo_lista, transistore_lista = zl3.elementu_ezlinealak(
        b, cir_el_extended)
    capacitor_lista, inductor_lista = zl4.elementu_dinamikoak(
        b, cir_el_extended)

    has_nonlinear = diodo_lista or transistore_lista
    has_dynamic = capacitor_lista or inductor_lista

    sorgailua, hasiera, amaiera, pausua, a, pr, op = zl2.tr_dc_parametroak(
        filename)

    if pr == ".pr":
        zl1.print_cir_info(cir_el_extended, cir_nd_extended, b, n,
                           nodes, el_num)
        print("\nIncidence Matrix: ")
        print(incidence_matrix)

    base_dir = os.path.dirname(filename)
    sims_dir = os.path.join(base_dir, "sims")
    os.makedirs(sims_dir, exist_ok=True)
    output_csv = os.path.join(sims_dir,
                              os.path.splitext(os.path.basename(filename))[0])


    if op == ".op":
        m = zl2.getMatrixes(b, cir_el_extended, cir_val_extended,
                            cir_ctr_extended)
        M, N, Us = m[0], m[1], m[2]

        if has_dynamic:
            Vc_prev = np.zeros(b)
            Il_prev = np.zeros(b)
            M, N, Us, Vc_prev, Il_prev = zl4.dynamic_NR(
                b, cir_el_extended, cir_val_extended, pausua, Vc_prev, Il_prev,
                capacitor_lista, inductor_lista, M, N, Us
            )

        if has_nonlinear:
            sol = zl3.Newton_Raphson(M, N, Us, A, At, n, b, cir_val_extended,
                                     diodo_lista, transistore_lista)
        else:
            T = zl2.get_T(A, At, M, N, Us, b, n)
            U = zl2.get_U(A, b, n, Us)
            zl2.SingularMatrix(T, U)
            sol = np.linalg.solve(T, U)

        zl2.print_solution(sol, b, n)

    if a == ".tr" or a == ".dc":
        if sorgailua != "0":
            output_csv += f"_{sorgailua}{a}"
        else:
            output_csv += a
        if has_dynamic:
            zl4.save_as_csv_p4(
                b, n, output_csv, cir_el_extended, cir_val_extended,
                cir_ctr_extended, A, At, cir_nd, cir_nd_extended,
                incidence_matrix, nodes, a, hasiera, amaiera, pausua, sorgailua
            )
        elif has_nonlinear:
            zl3.save_as_csv_p3(
                b, n, output_csv, cir_el_extended, cir_val_extended,
                cir_ctr_extended, A, At, cir_nd, cir_nd_extended,
                incidence_matrix, nodes, a, hasiera, amaiera, pausua, sorgailua
            )
        else:
            zl2.save_as_csv(
                b, n, output_csv, cir_el_extended, cir_val_extended,
                cir_ctr_extended, A, At, cir_nd, cir_nd_extended,
                incidence_matrix, nodes, a, hasiera, amaiera, pausua, sorgailua
            )
