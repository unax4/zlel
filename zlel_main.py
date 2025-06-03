#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: zlel_main.py
    :synopsis: Main module for circuit analysis with automatic
               analysis selection

.. moduleauthor:: Aitor Lavin (aitorlavin02@gmail.com)
                  Unax Arregi (arregiunax@gmail.com)
"""

import zlel.zlel_p1 as zl1
import zlel.zlel_p2 as zl2
import zlel.zlel_p3 as zl3
import zlel.zlel_p4 as zl4
import numpy as np
import time
import sys
import os

if __name__ == "__main__":
    start = time.perf_counter()

    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "Aitor_Lavin.Unax_Arregi/cirs/all/1_zlel_adibide_op.cir"

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

        end = time.perf_counter()
        # print(f"Elapsed time: {end - start} seconds")
