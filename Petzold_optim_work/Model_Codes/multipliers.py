# -*- coding: utf-8 -*-
"""
Multiplier vectors for ICs and parameters - set to ones for nominal realisation
Load vectors from module where needed (parameters.py, ICs.py)
"""

# V_c: vector for continuous parameters
# V_d: vector for discrete parameters (either 0 or 1, e.g. switches)
# V_IC: vector for initial conditions

# The parameters/ICs are multiplied by the elements in the vector

import numpy as np
import pandas as pd

# Nominal realisation (parameters and ICs just mulitplied by ones)
V_c = np.ones(234)
V_d = np.ones(72)
# V_d[0] not used
V_d[16]=1
V_IC = np.ones(51)

#%% Choose optimisation type
# optimised_once: the top 15 parameters are scaled by a multiplier (obtained by optimising the base model)
# optimised_twice_top15 = the top 15 parameters are scaled by a multiplier (obtained by optimising the 'optimised once' model)
# optimised_twice_all = all parameters are scaled by a mulitplier (obtained by optimising the 'optimised once' model)
# no_optimisation = nominal realisation

optimisation_type = 'optimised_twice_top25'


if optimisation_type  == 'optimised_once':
    # Optimised once: (from table in paper and original code)
    V_c[12] = 1.00522104739110 # k_syn
    V_c[52] = 1.00954355474713 # G_Na_k
    V_c[57] = 1.00398205815597 # J_NaK_max
    V_c[62] = 1.01824590459778 # v_6
    V_c[65] = 1.00757750864527 # R_s
    V_c[66] = 1.00065844660895 # R_k
    V_c[96] = 1.00225588624749 # K_Na_k
    V_c[99] = 1.01291065010860 # K_act
    V_c[102] = 1.00000325457168 # k_pump
    V_c[106] = 1.00383963867958 # K_G
    V_c[128] = 0.999991448721558 # v_Ca2_i
    V_c[129] = 1.00723565879070 # R_Ca_i
    V_c[178] = 1.00735943104210 # z_2
    V_c[179] = 1.00971655606907 # z_3
    V_c[181] = 1.00670962636797 # z_5
    
elif optimisation_type  == 'optimised_twice_top25':
    # Optimised twice: (from spreadsheet)
    V_c[12] = 0.901659341 # k_syn
    V_c[52] = 1.14190286 # G_Na_k
    V_c[57] = 1.075425867 # J_NaK_max
    V_c[60] = 1.099567279    # v_4
    V_c[62] = 1.093350074 # v_6
    V_c[65] = 0.9016100429 # R_s
    V_c[66] = 1.0987147 # R_k
    V_c[96] = 0.9717258082 # K_Na_k
    V_c[99] = 0.913333447 # K_act
    V_c[102] = 1.098161882 # k_pump
    V_c[106] = 1.09006668 # K_G
    V_c[117] = 1.09230102 # B_i
    V_c[125] = 1.026754194 # L_i
    V_c[126] = 1.011946741 # G_Ca_i
    V_c[127] = 0.974762433 # v_Ca1_i
    V_c[128] = 0.9833894888 # v_Ca2_i
    V_c[129] = 0.8429664916 # R_Ca_i
    V_c[136] = 1.099672269 #sigma_0
    V_c[141] = 0.9364855138 # G_k_i
    V_c[142] = 0.9684384341 # v_K_i
    V_c[175] = 1.002162089 # v_Ca3_i
    V_c[176] = 1.097160303 # R_K_i
    V_c[178] = 1.084332631 # z_2
    V_c[179] = 0.8927456823 # z_3
    V_c[181] = 1.010882651 # z_5
    
elif optimisation_type  == 'optimised_twice_all':
    # Import the multipliers for optimisation from the excel file (optimised twice, all parameters)
    data = pd.read_excel('./September2020ParameterNotes.xlsx')
    df = pd.DataFrame(data, columns = ['Final Optimization Modifier'])
    V_c = df.values.flatten()


