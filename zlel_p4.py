#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: zlel_p4.py
    :synopsis: Módulo para análisis de circuitos dinámicos con Euler Backward

.. moduleauthor:: Aitor Lavin (aitorlavin02@gmail.com)
                  Unax Arregi (arregiunax@gmail.com)
"""

import numpy as np
import sys
import math

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
    Calcula los parámetros del equivalente NR de un capacitor usando Euler Backward.
    
    Args:
        C: Capacitancia del condensador
        Vc_prev: Tensión en el condensador en el paso anterior
        h: Paso de tiempo
        
    Returns:
        gc: Conductancia equivalente
        Ic: Fuente de corriente equivalente
    """
    gc = 1 / (h / C)
    Ic = gc * Vc_prev
    return gc, Ic


def inductor_NR(L, Il_prev, h):
    """
    Calcula los parámetros del equivalente NR de un inductor usando Euler Backward.
    
    Args:
        L: Inductancia del inductor
        Il_prev: Corriente en el inductor en el paso anterior
        h: Paso de tiempo
        
    Returns:
        gl: Conductancia equivalente
        Il: Fuente de corriente equivalente
    """
    gl = h / L
    Vl = Il_prev * gl
    return gl, Vl


def elementu_dinamikoak(b, cir_el_extended):
    """
    Identifica elementos dinámicos (condensadores e inductores) en el circuito.
    
    Args:
        b: Número de ramas
        cir_el_extended: Array extendido de elementos del circuito
        
    Returns:
        Lista con dos sublistas: [condensadores, inductores]
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
    Updates M, N, Us matrices for dynamic elements using Backward Euler.
    """
    # Initialize current state variables
    Vc_actual = np.zeros(b)
    Il_actual = np.zeros(b)
    
    # Process capacitors
    for j in capacitor_lista:
        C = float(cir_val_extended[j][0])
        Vc0 = float(cir_val_extended[j][1]) if len(cir_val_extended[j]) > 1 else 0.0
        
        # Apply initial condition if first step
        if np.all(Vc_prev == 0) and Vc0 != 0:
            Vc_prev[j] = Vc0
            
        gc, Ic = capacitor_NR(C, Vc_prev[j], h)
        M[j,j] = gc
        N[j,j] = -1
        Us[j] = Ic
        Vc_actual[j] = Vc_prev[j]

    # Process inductors - corrected implementation
    for j in inductor_lista:
        L = float(cir_val_extended[j][0])
        Il0 = float(cir_val_extended[j][1]) if len(cir_val_extended[j]) > 1 else 0.0
        
        # Apply initial condition if first step
        if Il_prev[j] == 0 and Il0 != 0:
            Il_prev[j] = Il0
            
        gl, Vl = inductor_NR(L, Il_prev[j], h)
        M[j,j] = 1
        N[j,j] = -gl
        Us[j,0] = Vl  # Ensure Us is treated as 2D array
        Il_actual[j] = Il_prev[j]

    return M, N, Us, Vc_actual, Il_actual

def save_as_csv_p4(b, n, filename, cir_el_extended, cir_val_extended,
                  cir_ctr_extended, A, At, cir_nd, cir_nd_extended,
                  incidence_matrix, nodes, a, hasiera, amaiera, pausua,
                  sorgailua):
    """ 
    Transient simulation with dynamic elements.
    """
    header = zl2.build_csv_header("t", b, n)
    
    # Identify dynamic and nonlinear elements
    capacitor_lista, inductor_lista = elementu_dinamikoak(b, cir_el_extended)
    diodo_lista, transistore_lista = zl3.elementu_ezlinealak(b, cir_el_extended)
    
    # Initialize state variables
    Vc_prev = np.zeros(b)
    Il_prev = np.zeros(b)
    
    with open(filename, 'w') as file:
        print(header, file=file)
        
        t = hasiera
        while t < amaiera:
            # Get basic matrices
            if a == ".tr":
                m = zl2.op_tr(b, cir_el_extended, cir_val_extended,
                              cir_ctr_extended, t)
            elif a == ".dc":
                m = zl2.op_dc(b, cir_el_extended, cir_val_extended,
                              cir_ctr_extended, t, sorgailua)
            
            M, N, Us = m[0], m[1], m[2]
            
            # Check for errors
            zl1.erroreak(cir_nd_extended, cir_el_extended, cir_ctr_extended,
                         cir_val_extended, incidence_matrix, nodes, b)
            
            # Process dynamic elements
            M, N, Us, Vc_prev, Il_prev = dynamic_NR(
                b, cir_el_extended, cir_val_extended, pausua, 
                Vc_prev, Il_prev, capacitor_lista, inductor_lista, M, N, Us)
            
            # Solve nonlinear system
            sol = zl3.Newton_Raphson(M, N, Us, A, At, n, b, cir_val_extended,
                                    diodo_lista, transistore_lista)
            
            # Update state variables from solution
            for j in capacitor_lista:
                Vc_prev[j] = float(sol[n-1 + j])  # Explicit conversion to float
            
            for j in inductor_lista:
                # Correct inductor current update:
                # i_L = (v_L * h/L) + i_L_prev
                v_L = float(sol[n-1 + j])  # Voltage across inductor
                Il_prev[j] = (v_L * pausua / float(cir_val_extended[j][0])) + Il_prev[j]
            
            # Save results
            sol = np.insert(sol, 0, t)
            sol_csv = ','.join(['%.9f' % num for num in sol])
            print(sol_csv, file=file)
            
            t += pausua

if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "cirs/all/3_zlel_RC.cir"
    
    # Parsear el circuito
    cir_el, cir_nd, cir_val, cir_ctr = zl1.cir_parser(filename)
    modified_results = zl1.modify_function(cir_el, cir_nd, cir_val, cir_ctr)
    cir_el_extended = np.array(modified_results[0])
    cir_nd_extended = np.array(modified_results[1])
    cir_ctr_extended = np.array(modified_results[3])
    cir_val_extended = np.array(modified_results[2])
    
    # Obtener información básica del circuito
    nodes = zl1.get_nodes(cir_nd)
    b = zl1.get_branches(cir_el)
    n = zl1.get_number_of_nodes(nodes)
    el_num = zl1.get_elements(cir_el)
    
    # Construir matrices
    incidence_matrix = zl1.get_incidence_matrix(n, b, nodes, cir_nd_extended)
    A = zl1.get_incidence_murriztua(incidence_matrix, n)
    At = zl2.get_iraulia(A)
    
    # Identificar elementos dinámicos
    capacitor_lista, inductor_lista = elementu_dinamikoak(b, cir_el_extended)
    
    # Obtener parámetros de simulación
    sorgailua, hasiera, amaiera, pausua, a , pr, op= zl2.tr_dc_parametroak(filename)
    
    if pr==".pr":
        zl1.print_cir_info(cir_el_extended, cir_nd_extended, b, n, nodes, el_num)
        print("\nIncidence Matrix: ")
        print(zl1.get_incidence_matrix(n, b, nodes, cir_nd_extended))
    
    # Ejecutar simulación si es transitoria o DC
    if a == ".tr" or a == ".dc":
        output_csv = "3_zlel_arteztailea.csv"
        save_as_csv_p4(
            b, n, output_csv, cir_el_extended, cir_val_extended,
            cir_ctr_extended, A, At, cir_nd, cir_nd_extended,
            incidence_matrix, nodes, a, hasiera, amaiera, pausua, sorgailua
        )
        
        # Graficar resultados (ejemplo: tensión en nodo 1)
        zl2.plot_from_cvs(output_csv, "t", "e2", "Voltaje en nodo 2")
