# -*- coding: utf-8 -*-
"""
Initial conditions
"""

# Returns u0, an array containing the initial conditions for the state variables of the NVU model
# idx: dictionary containing the state variable names ('keys') and their corresponding index ('values')
# V_IC: vector used to modify the initial conditions (set to ones for normal simulation)
def set_initial_conditions(idx, V_IC):
    
    u0 = [None]*len(idx) # Initialise with correct size
    
    # Neuron
    u0[idx['E_t']] = 0 * V_IC[idx['E_t']]       # Proportion of excitatory cells firing per unit time [-] 
    u0[idx['I_t']] = 0 * V_IC[idx['I_t']]          # Proportion of inhibitory cells firing per unit time [-]
    u0[idx['K_e']] = 3.5 * V_IC[idx['K_e']]        # Extracellular K+ concentration [mM]
    u0[idx['Na_sa']] = 9.37 * V_IC[idx['Na_sa']]     # Na+ concentration in the soma/axon [mM]
    u0[idx['Na_d']] = 9.42 * V_IC[idx['Na_d']]      # Na+ concentration in the dendrite [mM]
    u0[idx['O2']] = 0.02566 * V_IC[idx['O2']]     # Tissue oxygen concentration [mM] 
    u0[idx['CBV']] = 1.204  * V_IC[idx['CBV']] #0.794      # Nondimensional cerebral blood volume [-] 
    u0[idx['HbR']] = 0.7641 * V_IC[idx['HbR']]     # Normalised deoxyhemoglobin concentration [-] 
    u0[idx['Ca_n']] = 0.1 * V_IC[idx['Ca_n']]       # Neuronal Ca2+ concentration [uM]
    u0[idx['nNOS']] = 0.01056 * V_IC[idx['nNOS']]   # Neuronal nitric oxide synthase concentration [uM]
    u0[idx['NO_n']] = 0.02425 * V_IC[idx['NO_n']]   # Neuronal nitric oxide concentration [uM] 
    
    # Astrocyte
    u0[idx['Na_k']] = 18740 * V_IC[idx['Na_k']]     # Astrocytic Na+ concentration [uM]
    u0[idx['K_k']] = 92660 * V_IC[idx['K_k']]      # Astrocytic K+ concentration [uM]  
    u0[idx['HCO3_k']] = 9085 * V_IC[idx['HCO3_k']]    # Astrocytic HCO3- concentration [uM]
    u0[idx['Cl_k']] = 8212 * V_IC[idx['Cl_k']]      # Astrocytic Cl- concentration [uM] 
    u0[idx['Na_s']] = 149200 * V_IC[idx['E_t']]    # Synaptic Na+ concentration [uM] 
    u0[idx['K_s']] = 2932 * V_IC[idx['K_s']]       # Synaptic K+ concentration [uM] 
    u0[idx['HCO3_s']] = 16980 * V_IC[idx['HCO3_s']]   # Synaptic HCO3- concentration [uM] 
    u0[idx['K_p']] = 3039 * V_IC[idx['K_p']]       # Perivascular K+ concentration [uM] 
    u0[idx['w_k']] = 8.26e-5 * V_IC[idx['w_k']]    # Open probability of the astrocytic BK channel [-] 
    u0[idx['Ca_k']] = 0.1435 * V_IC[idx['Ca_k']]    # Astrocytic Ca2+ concentration [uM] 
    u0[idx['s_k']] = 480.8 * V_IC[idx['s_k']]      # Ca2+ concentration in the astrocytic ER [uM] 
    u0[idx['h_k']] = 0.4107 * V_IC[idx['h_k']]     # Inactivation variable of the astrocytic IP3R channel [-] 
    u0[idx['I_k']] = 0.048299 * V_IC[idx['I_k']]   # Astrocytic IP3 concentration [uM]
    u0[idx['eet_k']] = 0.4350 * V_IC[idx['eet_k']]   # Astrocytic EET concentration [uM] 
    u0[idx['m_k']] = 0.513 * V_IC[idx['m_k']]      # Open probability of the astrocytic TRPV4 channel [-] 
    u0[idx['Ca_p']] = 1853 * V_IC[idx['Ca_p']]      # Perivascular Ca2+ concentration [uM] 
    u0[idx['NO_k']] = 0.02234 * V_IC[idx['NO_k']]   # Astrocytic nitric oxide concentration [uM] 
    u0[idx['v_k']] = -88.79 * V_IC[idx['v_k']]     # Astrocytic membrane potential [mV] 
    u0[idx['AA_k']] = 9.3 * V_IC[idx['AA_k']]       # Astrocytic arachidonic acid concentration [uM]
    
    # SMC
    u0[idx['Ca_i']] = 0.2641 * V_IC[idx['Ca_i']]    # SMC Ca2+ concentration [uM]
    u0[idx['s_i']] = 1.1686 * V_IC[idx['s_i']]     # Ca2+ concentration in the SR of the SMC [uM]
    u0[idx['v_i']] = -34.7 * V_IC[idx['v_i']]      # SMC membrane potential [mV]
    u0[idx['w_i']] = 0.2206 * V_IC[idx['w_i']]     # Open probability of the SMC BK channel [-]
    u0[idx['I_i']] = 0.275 * V_IC[idx['I_i']]      # SMC IP3 concentration [uM]
    u0[idx['NO_i']] = 0.02047 * V_IC[idx['NO_i']]   # SMC nitric oxide concentration [uM] 
    u0[idx['E_b']] = 0.6372 * V_IC[idx['E_b']]     # Fraction of sGC in the basal state [-]
    u0[idx['E_6c']] = 0.2606 * V_IC[idx['E_6c']]    # Fraction of sGC in the intermediate form [-] 
    u0[idx['cGMP_i']] = 6.1 * V_IC[idx['cGMP_i']]     # SMC cGMP concentration [uM]
    u0[idx['H_i']] = 0.096 * V_IC[idx['H_i']] #0.069   # SMC 20-HETE concentration [uM] 0.069
    u0[idx['AA_i']] = 9.3 * V_IC[idx['AA_i']]       # SMC arachidonic acid concentration [uM]
    
    # EC
    u0[idx['Ca_j']] = 0.8339 * V_IC[idx['Ca_j']]    # EC Ca2+ concentration [uM]
    u0[idx['s_j']] = 0.6262 * V_IC[idx['s_j']]     # Ca2+ concentration in the ER of the EC [uM]
    u0[idx['v_j']] = -68.39 * V_IC[idx['v_j']]     # EC membrane potential [mV]
    u0[idx['I_j']] = 0.825 * V_IC[idx['I_j']]      # EC IP3 concentration [uM]
    u0[idx['eNOS']] = 0.4451 * V_IC[idx['eNOS']]    # Endothelial nitric oxide synthase concentration [uM]
    u0[idx['NO_j']] = 0.02051 * V_IC[idx['NO_j']]   # EC nitric oxide concentration [uM] 
    
    # Wall Mechanics
    u0[idx['Mp']] = 0.0842 * V_IC[idx['Mp']]      # Fraction of free phosphorylated cross-bridges [-]
    u0[idx['AMp']] = 0.0622 * V_IC[idx['AMp']]     # Fraction of attached phosphorylated cross-bridges [-]
    u0[idx['AM']] = 0.2746 * V_IC[idx['AM']]      # Fraction of attached dephosphorylated cross-bridges [-]
    u0[idx['R']] = 22.44 * V_IC[idx['R']]        # Arteriolar radius [um]
    
    # GABA, Glu and NPY
#    u0[idx['GABA']] = 0         # Normalised GABA concentration [-]
#    u0[idx['Glu']] = 0          # Normalised glutamate concentration [-]
#    u0[idx['NPY']] = 0          # Normalised NPY concentration [-]
    return u0