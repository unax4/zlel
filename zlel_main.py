#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: zlel_main.py
    :synopsis:

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

"""
https://stackoverflow.com/questions/419163/what-does-if-name-main-do
https://stackoverflow.com/questions/19747371/
python-exit-commands-why-so-many-and-when-should-each-be-used
"""
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
    iraulia = zl2.get_iraulia(incidence_matrix)

    "INTZIDENTZI MATRIZE MURRIZTUA"

    A = zl1.get_incidence_murriztua(incidence_matrix, n)

    "INTZIDENTZI MATRIZE MURRIZTUAREN IRAULIA"

    At = zl2.get_iraulia(A)

    "PRAKTIKA HONETAKO FUNTZIOEN BALIOAK LORTU"

    m = zl2.getMatrixes(b, cir_el_extended, cir_val_extended, cir_ctr_extended)
    M = m[0]
    N = m[1]
    Us = m[2]
    izena=os.path.basename(filename)
    if izena[0]=="0":
        zl1.print_cir_info(cir_el_extended, cir_nd_extended, b, n, nodes, el_num)
        print("\nIncidence Matrix: ")
        print(incidence_matrix)
        zl1.erroreak(cir_nd_extended, cir_el_extended, cir_ctr_extended,
             cir_val_extended, incidence_matrix, nodes, b)

    elif izena[0]=="1":

        "w(T)"

        w = zl2.w_matrizea(nodes, b)

        "IDENTITATE MATRIZEA"

        identitate_matrizea = zl2.get_identitate_matrizea(M, At)

        "T ETA U"

        T = zl2.get_T(A, At, M, N, Us, b, n)
        U = zl2.get_U(A, b, n, Us)

        "ERROREAK AURKITU"

        zl1.erroreak(cir_nd_extended, cir_el_extended, cir_ctr_extended,
                    cir_val_extended, incidence_matrix, nodes, b)

        "PRAKTIKA HONETAKO ERROREA (SINGULAR MATRIX)"

        zl2.SingularMatrix(T, U)

        "SOLUZIOA"

        sol = np.linalg.solve(T, U)
        

        "sIMULAZIOAK"
        sorgailua, hasiera, amaiera, pausua, a , pr , op= zl2.tr_dc_parametroak(filename)
        if pr==".pr":
            zl1.print_cir_info(cir_el_extended, cir_nd_extended, b, n, nodes, el_num)
            print("\nIncidence Matrix: ")
            print(zl1.get_incidence_matrix(n, b, nodes, cir_nd_extended))

        if op==".op":
            zl2.print_solution(sol, b, n)

        if a == ".tr" or a == ".dc":
            # Ejecutar simulación y guardar resultados en CSV
            output_csv = izena+".csv"
            zl2.save_as_csv(
                b, n, output_csv, cir_el_extended, cir_val_extended,
                cir_ctr_extended, A, At, cir_nd, cir_nd_extended,
                incidence_matrix, nodes, a, hasiera, amaiera, pausua, sorgailua
                )

            # Graficar resultado de la simulación
            #zl2.plot_from_cvs(output_csv, "t", "e1", "Voltaje en nodo 1")

    elif izena[0]=="2":
        diodo_lista, transistore_lista = zl3.elementu_ezlinealak(b, cir_el_extended)
        sol = zl3.Newton_Raphson(M, N, Us, A, At, n, b, cir_val_extended,
            diodo_lista, transistore_lista)
        

        "Simulazioak"
        sorgailua, hasiera, amaiera, pausua, a , pr, op= zl2.tr_dc_parametroak(filename)
        if pr==".pr":
            zl1.print_cir_info(cir_el_extended, cir_nd_extended, b, n, nodes, el_num)
            print("\nIncidence Matrix: ")
            print(zl1.get_incidence_matrix(n, b, nodes, cir_nd_extended))

        if op==".op":
            zl2.print_solution(sol, b, n)

        if a == ".tr" or a == ".dc":
            # Ejecutar simulación y guardar resultados en CSV
            output_csv = izena+".csv"
            zl3.save_as_csv_p3(
                b, n, output_csv, cir_el_extended, cir_val_extended,
                cir_ctr_extended, A, At, cir_nd, cir_nd_extended,
                incidence_matrix, nodes, a, hasiera, amaiera, pausua, sorgailua
                )

            # Graficar resultado de la simulación
            #zl2.plot_from_cvs(output_csv, "t", "e1", "Voltaje en nodo 1")

    elif izena[0]=="3":
        sorgailua, hasiera, amaiera, pausua, a , pr, op= zl2.tr_dc_parametroak(filename)
        if pr==".pr":
            zl1.print_cir_info(cir_el_extended, cir_nd_extended, b, n, nodes, el_num)
            print("\nIncidence Matrix: ")
            print(zl1.get_incidence_matrix(n, b, nodes, cir_nd_extended))
        

        # Ejecutar simulación si es transitoria o DC
        if a == ".tr" or a == ".dc":
            output_csv = izena+".csv"
            zl4.save_as_csv_p4(
                b, n, output_csv, cir_el_extended, cir_val_extended,
                cir_ctr_extended, A, At, cir_nd, cir_nd_extended,
                incidence_matrix, nodes, a, hasiera, amaiera, pausua, sorgailua
            )
            
            # Graficar resultados (ejemplo: tensión en nodo 1)
        #zl2.plot_from_cvs(output_csv, "t", "e1", "Voltaje en nodo 1")
    else:
        print("Error: File name does not start with 1, 2 or 3.")
        sys.exit(1)