# -*- coding: utf-8 -*-
"""
Model parameters

This file can be loaded as a module (e.g. put "import parameters as p" then access parameters as p.startpulse, p.dt etc)
"""
class parameters():   
    pass

def p_function(iteration,switch,change_index):
    # The vectors that scale the parameters. V_c: continuous, V_d: discrete (i.e. either 0 or 1)
    from multipliers import V_c_function, V_d_function
    V_c=V_c_function(iteration)
    V_d=V_d_function(change_index)
    p = parameters()
    
    # Latest changes:   tau_MTT 3e3 -> 1e3
    #                   k_syn 10.5 -> 11.5
    #                   G_GABA .2 * G_Cl_i -> .24 * G_Cl_i
    #                   npy_increase 0.04 -> 0.06
    #                   Hshift 30 -> 10
    # Changed to fit with correctly scaled whisker/opto data (previously HbO and HbR were divided by 100, now HbO is divided by 70 and HbR by 30)
    
    ''' Different model cases'''
    p.InputCase = 'SquarePulse'            # 'SquarePulse': Model with P,Q as a square input pulse [default]
                                         # 'ThalamicSquarePulse': Thalamic input T (Pinto et al) as a square input pulse
                                         # 'Inhibitory': model with inhibitory neurons only (P = 0, Q as a square pulse)
                                         # 'ZhengData': Model with P,Q scaled by Zheng experimental data
                                         # '2SquarePulses': Model with 2 square input pulses following the Zheng paradigm (second pulse 2 sec long)
                                         # 'ZhengFittedParams': Model with P,Q scaled by Zheng experimental data and several parameters fitted to experiment (Yves)
                                         # 'ThalamicTriangles': Thalamic input T (Pinto et al) as a series of fast triangular pulses determined in T_pulses() (in model_functions module)
                                         # 'ThalamicTrianglesZheng': Thalamic input T (Pinto et al) as a series of fast triangular pulses determined in T_pulses() (in model_functions module) scaled by Zheng experimental data
    
    p.NeuronType = 'whisker'   # 'whisker': normal excitatory neuron (whisker stimulation, somatosensory cortex) [default]
                             # 'opto': inhibitory interneuron (optogenetic response), note that this isn't working properly for the Thalamic input cases yet!!
                             # 'whiskerandopto': combo of both whisker and opto
    
    p.NOswitch = switch     # 'normal': eNOS and nNOS production [default]
                            # '7NI': no nNOS production
                            # 'LNAME': no eNOS or nNOS production
                            
    p.HETswitch = 'normal'    # 'normal' 20_HETE production on
                            # 'HET' 20_HETE production off
           

                     
    if (p.NeuronType == 'whisker' or p.NeuronType == 'whiskerandopto') and p.InputCase != 'Inhibitory':
        p.Pinput = 1.0 
        p.Qinput = 0#1.0
    else: 
        p.Pinput = 0.0 # inhibitory input only, for NeuronType = 'opto' or special 'Inhibitory' whisker case
        p.Qinput = 1.0
                        
    if p.NOswitch == 'normal':
        p.NOswitch_NE = 1 * V_d[3]      # Turn on neuronal NO production 
        p.NOswitch_EC_WSS = 1 * V_d[4]  # Turn on WSS induced eNOS production 
        p.NOswitch_EC_CA = 1 * V_d[5]   # Turn on Ca2+ induced eNOS production 
    elif p.NOswitch == '7NI':
        p.NOswitch_NE = 0 * V_d[3]     # Turn off neuronal NO production 
        p.NOswitch_EC_WSS = 1 * V_d[4] # Turn on WSS induced eNOS production 
        p.NOswitch_EC_CA = 1 * V_d[5] # Turn on Ca2+ induced eNOS production 
    elif p.NOswitch == 'LNAME':
        p.NOswitch_NE = 0 * V_d[3]     # Turn off neuronal NO production 
        p.NOswitch_EC_WSS = 0 * V_d[4] # Turn off WSS induced eNOS production 
        p.NOswitch_EC_CA = 0 * V_d[5] # Turn off Ca2+ induced eNOS production
    
    if p.HETswitch == 'normal':
       p.HETswitch_20HETE = 1 # turn on 20-HETE production
    elif p.HETswitch == 'HET':
        p.HETswitch_20HETE =0 # turn off 20-HETE production
       
    ''' Input [ms] '''
    p.startpulse = 500e3       # Start time of input stimulation, set as 500 sec to allow the model to reach a steady state
    p.lengthpulse = 2e3        # Length of stimulation
    p.Tend = 550e3             # End of simulation
    p.dt = 0.001e3             # Time step (this value is overwritten to be small if InputCase = 'ThalamicTriangles' or 'ThalamicTrianglesZheng') [0.001e3 default]
    p.double_pulse = 0         # 1: use the second pulse in the Zheng input data, 0: only one pulse. Used for InputCase = 'ZhengData' or 'ZhengFittedParams' [0:default]
    
    ''' Neuron '''
    
    #if NeuronType == 'whisker': # less GABA and NPY as there is a smaller population of inhibitory neurons in the whisker case
    #    alpha_GABA = 2.6e-3     # used to scale inflow of GABA due to increase in I_t, M.E. (Maite)
    #    alpha_NPY = 2.6e-3      # used to scale inflow of NPY due to increase in I_t, M.E. (Maite)
    #elif NeuronType == 'opto' or NeuronType == 'whiskerandopto':
    #    alpha_GABA = 4.2e-3     # used to scale inflow of GABA due to increase in I_t, M.E. (Maite)
    #    alpha_NPY = 4.2e-3      # used to scale inflow of NPY due to increase in I_t, M.E. (Maite)    
    
    if p.InputCase == 'ThalamicSquarePulse' or p.InputCase == 'ThalamicTriangles' or p.InputCase == 'ThalamicTrianglesZheng': # Thalamic input (Pinto et al)
        p.EI_relative = 1.499         # [-]    Maximum value for abs(E-I) after external inputs (block pulse)
        p.EImin = 0.139               # [-]    Minimum value for abs(E-I)
        p.I_relative = 1.512          # [-]    Maximum value for I after external inputs (block pulse)
        p.Imin = 0.225                # [-]    Minimum value for I
    else: # P, Q input
        p.EI_relative = 0.2894#0.268         # [-]    Maximum value for abs(E-I) after external inputs
        p.EImin = 0                   # [-]    Minimum value for abs(E-I) 
        p.Imin = 0                    # [-]    Minimum value for I    
    #    if NeuronType == 'whisker' or NeuronType == 'whiskerandopto':
    #        I_relative = 0.179      # [-]    Maximum value for I after external inputs 
    #    elif NeuronType == 'opto':
    #        I_relative = 0.016      # [-]    Maximum value for I after external inputs, smaller maximum when only inhibitory neurons active (P=0)
    
    # Wilson and Cowan model parameters for K_e, Na_sa, Na_d
    # Alpha, beta and rest values adjusted for the different inputs
    if p.InputCase == 'SquarePulse' or p.InputCase == '2SquarePulses' or p.InputCase == 'Inhibitory':
        p.alpha_K_e = 2.1#2.0
        p.beta_K_e = 4.2e-3
        p.alpha_Na_sa = 4.23
        p.beta_Na_sa = 0.39e-3
        p.alpha_Na_d = -2.12
        p.beta_Na_d = 0.75e-3
        p.K_eBase = 3.5
        p.Na_saBase = 9.37
        p.Na_dBase = 9.42
        p.tau_e = 3.0      
        p.tau_i = 3.0           
    elif p.InputCase == 'ZhengData':
        p.alpha_K_e = 2.2
        p.beta_K_e = 4.2e-3
        p.alpha_Na_sa = 4.23
        p.beta_Na_sa = 0.65e-3
        p.alpha_Na_d = -2.12
        p.beta_Na_d = 0.7e-3
        p.K_eBase = 3.567
        p.Na_saBase = 9.237
        p.Na_dBase = 9.318
        p.tau_e = 3.0        
        p.tau_i = 3.0    
    elif p.InputCase == 'ThalamicSquarePulse' or p.InputCase == 'ThalamicTriangles' or p.InputCase == 'ThalamicTrianglesZheng': # Thalamic input
        p.alpha_K_e = 1.8
        p.beta_K_e = 4.2e-3 #1.7e-3
        p.alpha_Na_sa = 4.23
        p.beta_Na_sa = 0.39e-3
        p.alpha_Na_d = -2.12
        p.beta_Na_d = 0.75e-3
        p.K_eBase = 3.5
        p.Na_saBase = 9.37
        p.Na_dBase = 9.42
        p.tau_e = 5.0        
        p.tau_i = 15.0    
    elif p.InputCase == 'ZhengFittedParams':
        p.alpha_K_e = 3.3    # 3.3 with wallMech = 2.0 for CBF fit to Zheng, 3.0 with wallMech = 2.0 for fit to Berwick
        p.beta_K_e = 10e-3
        p.alpha_Na_sa = 4.23
        p.beta_Na_sa = 0.65e-3
        p.alpha_Na_d = -2.12
        p.beta_Na_d = 0.7e-3
        p.K_eBase = 3.567
        p.Na_saBase = 9.237
        p.Na_dBase = 9.318
        p.tau_e = 3        
        p.tau_i = 3  
    
     # Max size of thalamic input function T(t)
    if p.InputCase == 'ThalamicTrianglesZheng':   # multiple fast triangular pulses scaled by Zheng input data
        p.Tinput = 1
    elif p.InputCase == 'ThalamicTriangles':      # multiple fast triangular pulses
        p.Tinput = 0.8
    elif p.InputCase == 'ThalamicSquarePulse':    # square block pulse, effectively a time 'average' of a series of fast triangular pulses
        p.Tinput = 0.5
    else:
        p.Tinput = 0
        
    # Zheng 2010 paradigm
    if p.InputCase == 'ZhengData' or p.InputCase == 'ThalamicTrianglesZheng' or p.InputCase == '2SquarePulses' or p.InputCase == 'ZhengFittedParams':
        p.ISI = 8e3                                           # Inter stimulus interval (ISI), i.e. the time in ms between stimulations [8e3:default]
        p.ISI_array = [0.6e3, 1e3, 2e3, 3e3, 4e3, 6e3, 8e3]   # array of possible ISI
        p.ISI_idx = p.ISI_array.index(p.ISI)                      # INDEX for the ISI
        
        # stim = INDEX for length of initial stimulation [2,8,16] beginning at zero
        if p.lengthpulse == 2e3:
            p.stim = 0
        elif p.lengthpulse == 8e3:
            p.stim = 1
        elif p.lengthpulse == 16e3:
            p.stim = 2
        else:
            print("When using this InputCase you must have lengthpulse = 2, 8 or 16 sec")
    
    if p.InputCase == 'ZhengData' or p.InputCase == '2SquarePulses' or p.InputCase == 'ZhengFittedParams':
        p.secondpulse = p.startpulse + p.lengthpulse + p.ISI  
        p.secondlength = 2e3                       # length of second pulse
    else:
        p.ISI = 0       # No second pulse
        p.secondlength = 0
    
    if p.InputCase == 'ThalamicTriangles' or p.InputCase == 'ThalamicTrianglesZheng':
        p.dt = 0.001e3 # Thalamic input with fast triangles needs a small timestep
    
    # E and I parameters
    p.c1 = 12              # [-]
    p.c2 = 10              # [-]
    p.c3 = 13              # [-]
    p.c4 = 11              # [-]
    p.r_e = 1.0            # [ms]?
    p.r_i = 4.0            # [ms]?
    p.k_e_input = 10
    p.k_i_input = 10
    p.a_e = 1.2            # [-]
    p.theta_e = 2.8        # [-]
    p.a_i = 1.0            # [-]
    p.theta_i = 4.0        # [-]
    
    # Thalamic input parameters (Pinto et al)
    p.wee = 42 * V_c[3]
    p.wie = 24.6 * V_c[4]
    p.wei = 42 * V_c[5]
    p.wii = 18 * V_c[6]
    p.wte = 53.43 * V_c[7]
    p.wti = 68.4 * V_c[8]
    
    # Glutamate parameters
    p.Glu_max = 1846 * V_c[9] #******
    p.Glu_slope = 0.1 * V_c[10]      # [mM] Slope of sigmoidal (M.E.)
    p.Ke_switch = 5 * V_c[11]      # [mM] Threshold past which glutamate vesicle is released (M.E.)                               
    
    # Elshin model constants
    p.Farad = 96.485       # [C / mmol]
    p.ph = 26.6995         # [RT/F]
    p.k_syn = 10.5 * V_c[12]         # original vlue 10.5 [-] Coupling scaling factor, values can range between 1.06 to 14.95 according to estimation from experimental data of Ventura and Harris . Original 11.5
    
    # BOLD constants
    p.d = 0.4 * V_c[15]              # [-]
    p.a_1 = 3.4 * V_c[16]            # [-]   
    p.a_2 = 1 * V_c[17]              # [-]
    p.V_0 = 0.03 * V_c[18]           # [-]
    p.tau_MTT = 3e3 * V_c[13] #1e3        # [ms]  
    p.tau_TAT = 5e3 * V_c[14]       # [ms] 
    
    # Oxygen model
    p.alpha_O2 =  0.05 * V_c[19]     # [-] Percentage of ATP production indepent of O2 0.05
    p.K_init_e = 2.9 * V_c[20]       # [mM]
    p.Na_init_sa = 10 * V_c[21]      # [mM]
    p.Na_init_d = 10 * V_c[22]       # [mM]
    
    p.CBF_init = 3.2e-2 * V_c[23]    # [-]3.2e-2 *****
    p.gamma_O2 = 0.1 * V_c[24]       # [-] 0.1 
    
    # NO pathway 
    p.m_c = 4 * V_c[25]              # [-] 
    p.K_mA = 650 * V_c[26] #.352           # [-] - fit to Santucci2008 650/1846
    p.K_mB = 2800 * V_c[27] #1.517          # [-] - fit to Santucci2008 2800/1846
    p.v_n = -40 * V_c[28]            # [mV]  the neuronal membrane potential , assumed to be approx constant in this model
    p.G_M = 0.46 * V_c[29]           # original value 46e9 !! [pS] = [kg**-1*m**-2*ms**3*A**2] --> converted to ms
    p.P_Ca_P_M = 3.6 * V_c[30]       # [-]  the relative conductance of the NMDA channel to Ca2+ compared to monovalent ions
    p.Ca_ex = 2e3 * V_c[31]          # [uM]  the external calcium concentration (in Comerford+David2008: 1.5 mM!)
    p.M = 1.3e5 * V_c[32]            # [uM]  the concentration of monovalent ions in the neuron
    p.n_NR2A = 0.63 * V_c[33]        # [-]  average number of NR2A NMDA receptors per synapse (Santucci2008)
    p.n_NR2B = 11 * V_c[34]          # [-]  average number of NR2B NMDA receptors per synapse (Santucci2008)
    p.V_max_NO_n = 4.22e-3 * V_c[35] # [ms**-1]  maximum catalytic rate of NO production (Chen2006) - obtained from fig 6 & equ 17 & 18
    p.O2_n = 200 * V_c[36]           # [uM]  tissue O2 concentration in the neuron (M.E.)
    p.K_mO2_n = 243 * V_c[37]        # [uM]  Chen2006
    p.LArg_n = 100 * V_c[38]         # [uM]  
    p.K_mArg_n = 1.5 * V_c[39]       # [uM]      
    p.k_O2_n = 9.6e-9 * V_c[40]      # [uM**-2 ms**-1]  # (Kavdia2002)
    p.x_nk = 25 * V_c[41]            # [um]   (M.E.)
    p.D_cNO = 3300e-3 * V_c[42]      # [um**2 ms**-1]  Diffusion coefficient NO (Malinski1993)   
    p.V_spine = 8e-5 * V_c[43]       # [pL]  volume of the neuronal dritic spine Santucci2008
    p.k_ex = 1600e-3 * V_c[44]       # [ms**-1]  decay rate constant of internal calcium concentration Santucci2008
    p.Ca_rest = 0.1 * V_c[45]        # [uM]  resting calcium concentration (in Comerford+David2008: 2.830 mM in Santucci2008P: 0.1 \muM)
    p.lambda_buf = 20 * V_c[46]      # [-]  buffer capacity Santucci2008
    p.V_maxNOS = 0.697e-3 * V_c[47]     # [uM ms**-1]  M.E. original values 25e-6
    p.K_actNOS = 0.8 * V_c[48]   # [uM]  original values 9.27e-5
    p.mu2_n = 0.2e-2 * V_c[49]    # [ms**-1]  original values 0.0167e-3 
    
    # Switches
    p.GluSwitch = 1 * V_d[1]      # Turn on/off glutamate input [1:default]
    p.O2switch = 1 * V_d[2]        # O2switch = 0 ATP is plentiful, O2switch = 1 ATP is limited (oxygen-limited regime) [1:default]
    p.LCpathway = 0 + (1- V_d[8])        # Locus coeruleus pathway used for InputCase = 'ZhengData', 'ThalamicTrianglesZheng', 'ZhengFittedParams' [0:default]
    
    ''' Astrocyte '''
        
    # Conductances in uM mV**-1 ms**-1
    p.G_BK_k = 10.25e-3 * V_c[50]       
    p.G_K_k = 6907.77e-3 * V_c[51]    
    p.G_Na_k = 227.00e-3 * V_c[52] #original value 227
    p.G_NBC_k = 130.74e-3 * V_c[53]  
    p.G_KCC1_k = 1.728e-3 * V_c[54]
    p.G_NKCC1_k = 9.568e-3 * V_c[55]
    p.G_Cl_k = 151.93e-3 * V_c[56]      
    
    p.J_NaK_max = 2.38e1 * V_c[57]                        # original value 2.37e01 [uM ms**-1]    
    
    # Switches to turn on and off some things
    p.rhoSwitch = 1 * V_d[9] 
    p.trpv_switch = 1 * V_d[10] 
    
    p.rho_min = 0.1 * V_c[58]      
    p.rho_max = 0.7 * V_c[59]   
    
    p.v_4 = 8.0 * V_c[60]  # original value 8mV
    p.v_5 = 15 * V_c[61] #mV
    p.v_6 = -55 * V_c[62] #original value -55 mV 
    p.Ca_3 = 0.4 * V_c[63] # uM 
    p.Ca_4 = 0.35 * V_c[64] # uM 
    
    # Scaling Constants
    p.R_s = 2.79e-8 * V_c[65] # original value 2.79e-08 m
    p.R_k = 6.0e-8 * V_c[66] # original value 6.0e-08 m
    p.R_tot = 8.79e-8 * V_c[67] # m
    p.VR_sa = (p.R_s/p.R_k) * V_c[68]
    
    # Calcium in the Astrocyte Equations Constants
    p.delta = 1.235e-2 * V_c[69]
    p.VR_ER_cyt = 0.185 * V_c[70]
    
    p.J_max = 2880e-3 * V_c[98]                                  # [uM ms**-1] 
    p.K_I = 0.03 * V_c[79] # uM
    p.K_act = 0.17 * V_c[99] # original avlue 0.17uM
    p.k_on = 2e-3 * V_c[71] #8e-3                                      # [uM ms**-1]  previously 2e-3
    p.K_inh = 0.1 * V_c[72] #uM
    
    p.r_h = 4.8e-3 * V_c[73]                                     # [uM/ms] *****
    p.k_deg = 1.25e-3 * V_c[74]                                  # [ms**-1]
    p.V_eet = 72e-3 * V_c[75]                                    # [ms**-1]
    p.k_eet = 7.2e-3 * V_c[76]                                   # [ms**-1]
    p.Ca_k_min = 0.1 * V_c[77] # uM
    p.eet_shift = 2 * V_c[78] #mV/uM
    
    p.P_L = 0.0804e-3 * V_c[100]                                  # [uM ms**-1]
    p.V_max = 20e-3 * V_c[101]                                    # [uM ms**-1]
    p.k_pump = 0.24 * V_c[102] # original value 0.24uM
    
    #TRPV4
    p.Capmin_k = 2000 * V_c[80] #uM
    p.C_astr_k = 40 * V_c[81] #pF
    p.gam_cae_k = 200 * V_c[82] #uM
    p.gam_cai_k = 0.01 * V_c[83] #uM
    p.epshalf_k = 0.1 * V_c[84] # s
    p.kappa_k = 0.1 * V_c[85]
    p.v1_TRPV_k = 120 * V_c[86] #mV
    p.v2_TRPV_k = 13 * V_c[87] #mV
    p.t_TRPV_k = 0.9e3 * V_c[88]                                 # [ms] 
    p.Ca_decay_k = 0.5e-3 * V_c[89]                              # [ms**-1]
    p.G_TRPV_k = 3.15e-7 * V_c[90]                               # [ms**-1]
    p.r_buff = 0.05 * V_c[91] # Rate at which Ca2+ from the TRPV4 channel at the foot is buffered compared to rest of channels on the astrocyte body [-]
    
    # Perivascular space
    p.VR_pa = 0.001 * V_c[92] # [-]
    p.VR_ps = 0.001 * V_c[93] # [-]
    p.K_p_min = 3e3 * V_c[95] # uM
    p.R_decay = 0.15e-3 * V_c[94]                           # [ms**-1]  0.15e-3   ***0.15e-3*3
    
    # Fluxes Constants
    p.R_g = 8.315 #J mol**-1 K**-1 Gas constant
    p.T = 300 # K Temperature
    
    p.K_Na_k = 10000 * V_c[96] # original value 10000 uM
    p.K_K_s = 1500 * V_c[97]# uM
    
    # Additional Equations Astrocyte Constants
    p.z_K = 1# [-]
    p.z_Na = 1# [-]
    p.z_Cl = -1# [-]
    p.z_NBC = -1# [-]
    p.z_Ca = 2 # [-]
    p.BK_end = 40 * V_c[103] # [-]
    p.K_ex = 0.26 * V_c[104] #uM
    p.B_ex = 11.35 * V_c[105] #uM
    p.K_G = 8.82 * V_c[106] # original value 8,82
    p.psi_w = 2.664e-3 * V_c[107]                               # [ms**-1]
    
    # NO Pathway
    p.x_ki = 25 * V_c[108]            # [um]   (M.E.)
    p.k_O2_k = 9.6e-9 * V_c[109]      # [uM**-2 ms**-1]   (Kavdia2002)
    p.O2_k = 200 * V_c[110]           # [uM]   (M.E.)
    
    ''' SMC/EC '''
    
    # Smooth Muscle Cell ODE Constants
    p.gamma_i = 1970 #mV uM**-1
    p.lambda_i = 45e-3 * V_c[111]                                 # [ms**-1]
    
    # Endothelial Cell ODE Constants
    p.C_m_j = 25.8 * V_c[112]                                     # pF veranderen naar iets met ms??
    p.J_PLC = 0.11e-3 * V_c[113]                                 #0.11 for steady state, 0.3 for oscillations # [uM ms**-1] 
    p.J_0_j = 0.029e-3 * V_c[114]                                 # [uM ms**-1] constant Ca influx (EC)
    
    # Smooth Muscle Cell Flux Constants
    p.F_i = 0.23e-3 * V_c[115]                                    # [uM ms**-1]      # IP3/RYR channel strength
    p.K_r_i = 1 * V_c[116] #uM
    
    p.B_i = 2.08e-3 * V_c[117]                                   # original value 2.025e-3[uM ms**-1]     # SERCA pump strength     
    p.c_b_i = 1.0 * V_c[118] #uM
    
    p.C_i = 55e-3 * V_c[119]                                      # [uM ms**-1]
    p.s_c_i = 2.0 * V_c[120] #uM
    p.c_c_i = 0.9 * V_c[121] #uM
    
    p.D_i = 0.24e-3 * V_c[122]                                    # [ms**-1]
    p.v_d = -100 * V_c[123] #mV
    p.R_d_i = 250 * V_c[124] #mV
    
    p.L_i = 0.0256e-3 * V_c[125]                                   # original value 0.025e-3 ms**-1]
    
    p.G_Ca_i = 1.30e-6 * V_c[126]                                 # original value 1.29e-6 [uM mV**-1 ms**-1]
    p.v_Ca1_i = 100.0 * V_c[127] # original value 100mV
    p.v_Ca2_i = -24.0 * V_c[128]  #original avlue -24 mV
    p.R_Ca_i = 8.5 * V_c[129] # original value 8.5mV
    
    p.G_NaCa_i = 3.16e-6 * V_c[130]                              # [uM mV**-1 ms**-1]
    p.c_NaCa_i = 0.5 * V_c[131] #uM
    p.v_NaCa_i = -30 * V_c[132] #mV
    
    p.G_stretch = 6.1e-6 * V_c[133]                               # [uM mV**-1 ms**-1]   (Also EC parameter)
    p.alpha_stretch = 7.4e-3 * V_c[134] # mmHg**-1     (Also EC parameter)
    p.trans_p_mmHg = 30 * V_c[135] # mmHg             (Also EC parameter) transmural pressure. 30 mmHg = 4000 Pa
    p.sigma_0 = 500 * V_c[136] # original value 500 mmHg                 (Also EC parameter)
    p.E_SAC = -18 * V_c[137] # mV                     (Also EC parameter)
    
    p.F_NaK_i = 4.32e-5 * V_c[138]                               # [uM ms**-1]
    
    p.G_Cl_i = 1.34e-6 * V_c[139]                                 # [uM mV**-1 ms**-1]
    p.v_Cl_i = -25 * V_c[140] #mV
    
    p.G_K_i = 4.46e-6 * V_c[141]                                  # original value 4.46e-6[uM mV**-1 ms**-1]
    p.v_K_i = -94 * V_c[142] # original  value -94 mV
    
    p.F_KIR_i = 1.285e-9 * V_c[143] #0.381                       # [uM mV**-1 ms**-1] F_KIR_i = 1.285e-6
    p.k_d_i = 0.1e-3 * V_c[144]                                   # [ms**-1]
    
    # othelial Cell Flux Constants
    p.F_j = 0.23e-3 * V_c[145]                                    # [uM ms**-1]
    p.K_r_j = 1 * V_c[146] #uM
    p.B_j = 0.5e-3 * V_c[147]                                     # [uM ms**-1]
    p.c_b_j = 1 * V_c[148] #uM
    p.C_j = 5e-3 * V_c[149]                                       # [uM ms**-1]
    p.s_c_j = 2 * V_c[150] #uM
    p.c_c_j = 0.9 * V_c[151] #uM
    p.D_j = 0.24e-3 * V_c[152]                                    # [ms**-1]
    
    # (G_stretch, alpha_stretch, trans_p_mmHg, sigma0, E_SAC are included above in 
    #  SMC flux Constants)
    
    p.L_j = 0.025e-3 * V_c[153]                                   # [ms**-1] 
    
    p.G_cat_j = 6.6e-7 * V_c[154]                                 # [uM mV**-1 ms**-1]
    p.E_Ca_j = 50 * V_c[155] #mV
    p.m_3_cat_j = -0.18 * V_c[156] 
    p.m_4_cat_j = 0.37 * V_c[157] 
    
    p.G_tot_j = 6927e9 * V_c[158] #6927                               # [mpS] = [kg**-1*m**-2*ms**3*A**2] --> WRONG?
    p.v_K_j = -80 * V_c[159] #mV
    
    p.c = -0.4 * V_c[160]
    p.bb_j = -80.8 * V_c[161] #mV
    p.a_1_j = 53.3 * V_c[162] #mV
    p.a_2_j = 53.3 * V_c[163] # mV 
    p.m_3b_j = 1.32e-3 * V_c[164] #mV**-1
    p.m_4b_j = 0.3 * V_c[165] #mV
    p.m_3s_j = -0.28 * V_c[166]
    p.m_4s_j = 0.389 * V_c[167] 
    
    p.G_R_j = 955e9 * V_c[168] #e9                                    # [mpS] = [ ] --> WRONG?
    p.v_rest_j = -31.1 * V_c[169] #mV
    p.k_d_j = 0.1e-3 * V_c[170]                                   # [ms**-1]
    
    p.P_Ca = 0.05e-3 * V_c[171]                                   # [ms**-1]
    p.P_IP3 = 0.05e-3 * V_c[172]                                  # [ms**-1]
    p.G_coup = 0.5e-3 * V_c[173]                                  # [ms**-1]
    
    # Additional Equations Constants
    p.alpha_act_i = 0.13 * V_c[174] #uM**2
    p.v_Ca3_i = -27.0 * V_c[175] #original value -27mV
    p.R_K_i = 12.0 * V_c[176] # original value 12mV
    p.z_1 = 4.5e-3 * V_c[177] #mV uM**-1
    p.z_2 = 112 * V_c[178]  # original value 112mV
    p.z_3 = 4.2e-4 * V_c[179] # original value 4.2e-4 uM**-1
    p.z_4 = 12.6 * V_c[180] # -  
    p.z_5 = -7.4e-2 * V_c[181]  # original value -7.4e-2 mV**-1
    
    # NO pathway
    p.K_mArg_j = 1.5 * V_c[182] # [uM] 
    p.K_mO2_j = 7.7 * V_c[183] # [uM]  Chen2006
    p.k_dno = 0.01e-3 * V_c[184]                                  # [ms**-1] 
    p.K_m_mlcp = 5.5 * V_c[185] # [uM] 
    p.V_NOj_max = 1.22e-3 * V_c[186]                             # [ms**-1]  
    p.LArg_j = 100 * V_c[187] # [uM] 
    p.k_O2 = 9.6e-9 * V_c[188]                                    # [uM**-2 ms**-1] 
    p.W_0 = 1.4 * V_c[189] #  shear gating constant (Comerford2008)
    p.delta_wss = 2.86 * V_c[190] # [Pa]  the membrane shear modulus (Comerford2008)
    p.k_1 = 100e-3 * V_c[191]                                     # [ms**{-1}] 
    p.k1 = 2 * V_c[192]                                           # [uM**-1 ms**-1] 
    p.k2 = 0.1e-3 * V_c[193]                                      # [ms**-1] 
    p.k3 = 3e-3 * V_c[194]                                        # [uM**-1 ms**-1] 
    p.V_max_sGC = 0.8520e-3 * V_c[195]                            # [uM/ms] 
    p.k_pde = 0.0195e-3 * V_c[196]                                # [ms**-1] 
    p.C_4 = 0.011e-3 * V_c[197]                                   # [ms**{-1} microM**{-2}] 
    p.K_m_pde = 2 * V_c[198] # [uM] 
    p.gam_eNOS = 0.1 * V_c[199] # [-] 
    p.mu2_j = 0.0167e-3 * V_c[200]                                # [ms**-1] 
    p.K_dis = 9e-5 * V_c[201]                                     # [uM ms**-1] 
    p.K_eNOS = 4.5e-1 * V_c[202] # [uM]     
    p.g_max = 0.06e-3 * V_c[203]                                 # [uM ms**-1]     
    p.alp = 2 * V_c[204] # [-]  zero shear open channel constant (Comerford2008 in Wiesner1997: alp = 3
    p.delta_p_L = 9.1e-2 * V_c[205] # Pa/um ME
    p.x_ij = 3.75 * V_c[206] # [um]  (Kavdia2002)
    
    ''' Wall mechanics '''
    p.wallMech = 8.7 * V_c[0]                                  # ****scaling factor for the wall mechanics, original = 2.7 
    
    # Contraction Equation Constants
    p.K_3 = 0.4e-3 * V_c[207]                                     # [ms**-1]
    p.K_4 = 0.1e-3 * V_c[208]                                     # [ms**-1]
    p.K_7 = 0.1e-3 * V_c[209]                                     # [ms**-1]
    p.gamma_cross = 17e-3 * V_c[210]                              # [uM**-3 ms**-1]
    p.n_cross = 3 * V_c[211]                                      # [-]
    
    # Mechanical Equation Constants
    p.eta_R = 1e7 * V_c[212]                                      # [Pa ms]
    p.R_init = 20 * V_c[213]                                      # [um] 
    p.trans_p = 4000 * V_c[214]                                   # [Pa]
    p.E_passive = 66e3 * V_c[215]                                 # [Pa]
    p.E_active = 233e3 * V_c[216]                                 # [Pa]
    p.alpha = 0.6 * V_c[217]                                      # [-]
    p.k_mlcp_b = 0.0086e-3 * V_c[218]                             # [ms**-1]
    p.k_mlcp_c = 0.0327e-3 * V_c[219]                             # [ms**-1]
    
    ''' 20-HETE parameters ''' #*** we should probably adjust these so 20-HETE doesn't have such a huge effect
    
    p.AA_max = 29 * V_c[220]  # uM
    p.AA_m = 0.161 * V_c[221]  # uM
    p.Ca0 = 0.1432 * V_c[222]  # uM
    p.D_AA = .033 * V_c[223] #.5*.033    # um**2/ms, model estimate
    p.tau_AA = (p.x_ki ** 2 /  (2 * p.D_AA)) * V_c[224] # Time constant for AA diffusion between astrocyte and SMC
    
    p.V_a = 0.212e-2 * V_c[225]  # ms**-1  old value 0.212e-3
    p.K_a = 228.2 * V_c[226]  # uM
    p.V_f = 0.0319e-2 * V_c[227]  # ms**-1 old value 0.0319e-3
    p.K_f = 23.5 * V_c[228]  # uM
    p.lambda_h = 2.0e-3 * V_c[229]   # ms**-1  old value 0.139e-3
    p.Hshift = 30 * V_c[230]  #10
    p.H0 = 0.126 * V_c[231] #0.068   # 0.126 was wrong baseline value before, must have changed when other parameters changed
    
    p.NO_rest = 0.02047 * V_c[232]
    p.R_NO = 0.02 * V_c[233]
    
    ### old param
    p.CBF_0 = 5.403e-2
    ###
    
    
    #''' Oxygen extraction fraction parameters '''
    #
    #f_in0 = 1.521           # Baseline f_in for normalising, note that CBF/CBF_init is not exactly 1 so we need to do this to get f_in = 1 at baseline #1.681
    p.O2_0 = 0.01             # previously 0.02 (from Chang2013) [mM] baseline tissue O2 concentration
    p.g_OEF = 0.2             # Ratio of tissue O2 to plasma O2 (0.01 / 0.053)
    p.E_0 = 0.4               # [-] baseline OEF
    p.J_0 = 0.053             # previously 0.032 [mM/s] steady state change in pump & background flux due to CBF
    #
    #''' GABA parameters '''
    #GABAswitch = 1          # 1: GABA (opens chlorine channels on AC and SMC) and no neuronal activity
    #
    #G_Tmin = 1              # [-] Min value of G_Tact (GABA-T activity), occurs when NO concentration is normal
    #G_Tmax = 2              # [-] Max value of G_Tact (GABA-T activity), occurs when NO = 0
    #GT_midpoint = 0.007     # [uM] midpoint of NO dependent G_Tact sigmoidal, M.E.
    #GT_slope = 0.003        # [uM] slope of NO dependent G_Tact sigmoidal, M.E.
    #
    #g_slope = 0.15          # [-] Slope of GABA sigmoidal, M.E.
    #g_midpoint = 0.8        # [-] Midpoint of GABA sigmoidal, M.E.  
    #E_GABA = -75            # [mV] Reversal potential of GABA activated Cl channel
    #G_GABA = .24 * G_Cl_i    # Maximum conductance of GABA activated Cl channel, based on SMC Cl leak channel conductance, M.E.
    #
    #beta_Glu = 4.2e-3       # clearance of Glu, M.E. 
    #beta_GABA = 4.2e-3      # clearance of GABA, M.E. (Maite)
    #GABAbase = 0            # baseline GABA concentration (Maite)
    #
    #''' NPY parameters '''
    #NPYswitch = 1           # 1: NPY release (opens VOCCs and causes constriction)  
    #
    #npy_increase = 0.06     # [-] proportion increase of VOCC conductance from baseline value due to NPY, e.g. 0.06 means 6% increase, M.E. 
    #npy_slope = 0.15        # [-] Slope of NPY sigmoidal, M.E. 
    #npy_midpoint = 0.8      # [-] Midpoint of NPY sigmoidal, M.E.
    #
    #beta_NPY = 4.2e-3       # clearance of NPY, M.E. (Maite)
    #NPYbase = 0             # baseline NPY concentration (Maite)
    #
    '''New stuff'''
    ## CICR channel in the astrocyte
    p.C_k = 30e-3   # [uM ms**-1] 5e-3 (EC), 55e-3 (SMC)
    p.s_c_k = 2    #uM
    p.c_c_k = 0.9  #uM
    
    return p