# -*- coding: utf-8 -*-
"""
Import the matlab data files
"""

from scipy import io, interpolate, stats
import numpy as np
import parameters as p
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind_from_stats

# Define a class that will contain data (since there are a lot of different arrays in some cases)
class Data():   
    pass

'''Import the Zheng 2010 neural input data'''
# Returns input_data, an array used for InputCase = 'ZhengData', 'ThalamicTrianglesZheng' or 'ZhengFittedParams' for scaling the neuronal stimulus input functions P, Q or T. Based on experimental data from Zheng 2010
# t: time vector
def import_Zheng_neural_data(t):
 
    # If the neuronal stimulation input profile is found from experimental data:
    if p.InputCase == 'ZhengData' or p.InputCase == 'ThalamicTrianglesZheng' or p.InputCase == 'ZhengFittedParams':
        
        # Import the data from mat file
        mat = io.loadmat('./mat_files/Zheng2010_data.mat')
        
        # Extract and flatten data
        neural_tim_vector = mat['neural_tim_vector'].flatten()
        neural_data = mat['neural_data']
   
        sum_neural_wh = [0]*len(neural_tim_vector) # Initialise array
        for animal in range(11):
            for experiment in range(10):
                sum_neural_wh = sum_neural_wh + neural_data[:, p.ISI_idx, p.stim, experiment, animal] # Sum all animal/trial data (110 total)
        mean_neural_wh = sum_neural_wh/110  # Average the neural data
    
        neural_tim_vector_shifted = neural_tim_vector*1e3 + p.startpulse - 20    # Shift so stimulation begins at p.startpulse, -20 so initial spike isn't chopped off
    
        # Interpolate so there is data for all timesteps for NVU
        interp_neural_f = interpolate.interp1d(neural_tim_vector_shifted, mean_neural_wh, fill_value = 0, bounds_error=False)
        interp_neural_wh = interp_neural_f(t)  
    
        if p.double_pulse == 0:  # Remove the second pulse if not wanted
            bool_array = np.in1d(t, p.startpulse + p.lengthpulse+1e3) # find index for where the second pulse begins
            t_idx = np.where(bool_array == True)
            t_idx = t_idx[0][0]
            interp_neural_wh[t_idx:] = 0
    
        # Locus coeruleus pathway - increases input stimulus after some period of time [default = 0]
        if p.LCpathway == 1:
            time_1 = np.linspace(0, p.lengthpulse, 10000) 
            I_LC = 0.0004 * (time_1/1e3)**2 + 1            # Chosen for the shape! Can be modified if needed
            time_1 = time_1 + p.startpulse                 # Shift so starts at p.startpulse
            
            # Interpolate so there is data for all timesteps for NVU
            I_LC_f = interpolate.interp1d(time_1, I_LC, fill_value = 1, bounds_error=False)
            I_LC = I_LC_f(t)  
        else:
             I_LC = 1
        
        input_data = interp_neural_wh * I_LC   # Save input profile for use in P(t), Q(t) or T(t)
                  
    else: # Neuronal stimulation input profile is just a square wave
        input_data = [] # Input data not used
        
    return input_data


'''Import the CBF data from Zheng 2010 experiments'''
# Returns:
# fig - figure of the averaged CBF data
# cbf_time - time variable
# cbf_Zheng - CBF variable
def import_Zheng_cbf_data():
 
    # If the neuronal stimulation input profile is found from experimental data:
    if p.InputCase == 'ZhengData' or p.InputCase == 'ThalamicTrianglesZheng' or p.InputCase == 'ZhengFittedParams':
        
        # Import the data from mat file
        mat = io.loadmat('./mat_files/Zheng2010_data.mat')
        
        # Extract and flatten data
        cbf_tim_vector = mat['cbf_tim_vector'].flatten()
        cbf_data = mat['cbf_data']
       
        sum_cbf = [0]*len(cbf_tim_vector) # Initialise array
        for animal in range(11):
            for experiment in range(10):
                sum_cbf = sum_cbf + cbf_data[:, p.ISI_idx, p.stim, experiment, animal] # Sum all animal/trial data (110 total)
        cbf_Zheng = sum_cbf/110 - 1  # Average the neural data
    
        cbf_time = cbf_tim_vector + p.startpulse*1e-3    # Shift so stimulation begins at p.startpulse, -20 so initial spike isn't chopped off

        fig = plt.figure(figsize=(9,6), dpi=80, edgecolor='k') 
        plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
        plt.xlabel('Time [s]')   
        plt.title('CBF (Zheng data)')
        plt.plot(cbf_time,cbf_Zheng)  
       
    else:
        print("Can only import Zheng CBF data if the neuronal stimulation input profile is taken from Zheng data (p.InputCase = 'ZhengData', 'ThalamicTrianglesZheng' or 'ZhengFittedParams')")
        fig, cbf_time, cbf_Zheng = plt.figure(), 0, 0

    return fig, cbf_time, cbf_Zheng

'''Import the data from the Berwick experiments using both air and oxygen'''
# ORIGINALLY CALLED two_sec_stim.mat !!!!!!!!!!!
# Returns:
# fig_airOxyComparison - figure of the air and oxygen experiments
# d - class containing the time and hemodynamic variables
# area - choose from: 'All', 'Whisker', 'Artery', 'Vein', 'Parenchyma'. Default = 'All'
def import_Berwick_AirOxyData(area='All'):
    
    mat = io.loadmat('./mat_files/Berwick_AirOxy_Data_Full.mat')
    
    # 19 = is the number of sessions (this data is from 5 animals the first four have 4 sessions the last animal has 3 - this is all in order (Animal 1=1:4, Animal 2=5 to 8) - There is one month gap between each of the sessions.
    # 4 = regions from where the responses are from (1 = whisker region, 2 = Artery, 3 = Vein, 4 = Parenchyma)
    # 200 = data length (time)
    # 30 = the number of trials in each of the 19 sessions
       
    # Initialise class to contain the data
    d = Data()
    
    d.time_airOxy = mat['tim'].flatten() - 5                # time vector
    
    '''Find the mean and standard deviation and plot to figure'''
    
    regions = { 'Whisker': 0, 'Artery': 1, 'Vein': 2, 'Parenchyma': 3} # Possible regions
    
    HbOair_array = mat['Hbo_2s_air_Fract']                  # extract the data from the mat file into a 4D array         
    HbRair_array = mat['Hbr_2s_air_Fract']
    HbTair_array = mat['Hbt_2s_air_Fract']
    HbOoxy_array = mat['Hbo_2s_oxy_fract']
    HbRoxy_array = mat['Hbr_2s_oxy_fract']
    HbToxy_array = mat['Hbt_2s_oxy_fract']
    
    if area == 'All':
        HbOair_array = HbOair_array.transpose(3,1,0,2)          # rearrange data so that it can be reshaped in the next step
        HbOair_array = np.reshape(HbOair_array, (19*4*30, 200)) # reshape into a 2D array so that the mean and std can be calculated
        HbRair_array = HbRair_array.transpose(3,1,0,2)
        HbRair_array = np.reshape(HbRair_array, (19*4*30, 200))   
        HbTair_array = HbTair_array.transpose(3,1,0,2)
        HbTair_array = np.reshape(HbTair_array, (19*4*30, 200))   
        HbOoxy_array = HbOoxy_array.transpose(3,1,0,2)
        HbOoxy_array = np.reshape(HbOoxy_array, (19*4*30, 200))         
        HbRoxy_array = HbRoxy_array.transpose(3,1,0,2)
        HbRoxy_array = np.reshape(HbRoxy_array, (19*4*30, 200))   
        HbToxy_array = HbToxy_array.transpose(3,1,0,2)
        HbToxy_array = np.reshape(HbToxy_array, (19*4*30, 200))   
    else:
        HbOair_array = HbOair_array[:,regions[area],:,:].transpose(2,0,1)   # extract data for specific area then rearrange data so that it can be reshaped in the next step
        HbOair_array = np.reshape(HbOair_array, (19*30, 200))               # reshape into a 2D array so that the mean and std can be calculated
        HbRair_array = HbRair_array[:,regions[area],:,:].transpose(2,0,1)   
        HbRair_array = np.reshape(HbRair_array, (19*30, 200))              
        HbTair_array = HbTair_array[:,regions[area],:,:].transpose(2,0,1)   
        HbTair_array = np.reshape(HbTair_array, (19*30, 200))               
        HbOoxy_array = HbOoxy_array[:,regions[area],:,:].transpose(2,0,1)  
        HbOoxy_array = np.reshape(HbOoxy_array, (19*30, 200))               
        HbRoxy_array = HbRoxy_array[:,regions[area],:,:].transpose(2,0,1)   
        HbRoxy_array = np.reshape(HbRoxy_array, (19*30, 200))              
        HbToxy_array = HbToxy_array[:,regions[area],:,:].transpose(2,0,1)  
        HbToxy_array = np.reshape(HbToxy_array, (19*30, 200))               
    
    d.HbOair_mean = np.mean(HbOair_array, axis=0)           # find the mean Hb profile
    d.HbOair_std = np.std(HbOair_array, axis=0)             # find the standard deviation for all the Hb profiles
    d.HbRair_mean = np.mean(HbRair_array, axis=0)
    d.HbRair_std = np.std(HbRair_array, axis=0)
    d.HbTair_mean = np.mean(HbTair_array, axis=0)
    d.HbTair_std = np.std(HbTair_array, axis=0)
    d.HbOoxy_mean = np.mean(HbOoxy_array, axis=0)
    d.HbOoxy_std = np.std(HbOoxy_array, axis=0)
    d.HbRoxy_mean = np.mean(HbRoxy_array, axis=0)
    d.HbRoxy_std = np.std(HbRoxy_array, axis=0)
    d.HbToxy_mean = np.mean(HbToxy_array, axis=0)
    d.HbToxy_std = np.std(HbToxy_array, axis=0)
    
    # Find the standard error
    d.HbOair_stderr = stats.sem(HbOair_array)             
    d.HbRair_stderr = stats.sem(HbRair_array)
    d.HbTair_stderr = stats.sem(HbTair_array)
    d.HbOoxy_stderr = stats.sem(HbOoxy_array)
    d.HbRoxy_stderr = stats.sem(HbRoxy_array)
    d.HbToxy_stderr = stats.sem(HbToxy_array)
    
    fig_airOxyComparison = plt.figure(figsize=(12,5), dpi=80, edgecolor='k') 
    fig_airOxyComparison.suptitle("Berwick data ({0:s})".format(area))
    
    # Find ylims based on max/min error bar for either graph
#    ylim1 = min(np.minimum(d.HbRair_mean - d.HbRair_std, d.HbRoxy_mean - d.HbRoxy_std)) - 0.005
#    ylim2 = max(np.maximum(d.HbOair_mean + d.HbOair_std, d.HbOoxy_mean + d.HbOoxy_std)) + 0.005
    ylim1 = 0.91
    ylim2 = 1.12
    
    ax = fig_airOxyComparison.add_subplot(1,2,1)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Air')  
    plt.errorbar(d.time_airOxy, d.HbOair_mean, d.HbOair_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time_airOxy, d.HbRair_mean, d.HbRair_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time_airOxy, d.HbTair_mean, d.HbTair_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.ylim(ylim1, ylim2)
    plt.legend()
    
    ax = fig_airOxyComparison.add_subplot(1,2,2)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Oxygen')  
    plt.errorbar(d.time_airOxy, d.HbOoxy_mean, d.HbOoxy_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time_airOxy, d.HbRoxy_mean, d.HbRoxy_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time_airOxy, d.HbToxy_mean, d.HbToxy_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.ylim(ylim1, ylim2)
    plt.legend()
    
    return fig_airOxyComparison, d
      
'''Import the LNAME/7NI data from Berwick et al'''
# Returns:
# fig1 - figure of the whisker experiments
# fig2 - figure of the optogenetic experiments
# d - class containing the time and hemodynamic variables
# area - choose from: 'All', 'Whisker', 'Artery', 'Vein', 'Parenchyma'. Default = 'All'
def import_Berwick_LNAME_7NI_data(area='All'):
    
    mat = io.loadmat('./mat_files/LNAME_data_preliminary.mat')
    
    # 8 = is the number of animals.
    # 4 = regions from where the responses are from (1 = whisker region, 2 = Artery, 3 = Vein, 4 = Parenchyma)
    # 200 = data length (time)
    
    # Initialise class to contain the data
    d = Data()
    
    d.time_LNAME = mat['tim'].flatten() - 5  # Time array
    
    regions = { 'Whisker': 0, 'Artery': 1, 'Vein': 2, 'Parenchyma': 3} # Possible regions
    saturation = { 'Whisker': 70, 'Artery': 80, 'Vein': 60, 'Parenchyma': 70} # saturation for each region, used for scaling the micromolar changes into fractional changes

    # Initialise arrays to contain the data for the 4 regions
    HbO_wh_preLNAME_all = [0]*4
    HbO_wh_postLNAME_all = [0]*4
    HbO_op_preLNAME_all = [0]*4
    HbO_op_postLNAME_all = [0]*4
    HbR_wh_preLNAME_all = [0]*4
    HbR_wh_postLNAME_all = [0]*4
    HbR_op_preLNAME_all = [0]*4
    HbR_op_postLNAME_all = [0]*4 
    HbT_wh_preLNAME_all = [0]*4
    HbT_wh_postLNAME_all = [0]*4
    HbT_op_preLNAME_all = [0]*4
    HbT_op_postLNAME_all = [0]*4
    
    # Extract and scale data for each region and each experiment type
    for r in regions:
        HbO_wh_preLNAME_all[regions[r]] = mat['hbo_whisker_pre'][:,regions[r],:]/saturation[r] + 1 
        HbO_wh_postLNAME_all[regions[r]] = mat['hbo_whisker_post'][:,regions[r],:]/saturation[r] + 1   
        HbO_op_preLNAME_all[regions[r]] = mat['hbo_opto_pre'][:,regions[r],:]/saturation[r] + 1   
        HbO_op_postLNAME_all[regions[r]] = mat['hbo_opto_post'][:,regions[r],:]/saturation[r] + 1     
        HbR_wh_preLNAME_all[regions[r]] = mat['hbr_whisker_pre'][:,regions[r],:]/(100-saturation[r]) + 1   
        HbR_wh_postLNAME_all[regions[r]] = mat['hbr_whisker_post'][:,regions[r],:]/(100-saturation[r]) + 1       
        HbR_op_preLNAME_all[regions[r]] = mat['hbr_opto_pre'][:,regions[r],:]/(100-saturation[r]) + 1   
        HbR_op_postLNAME_all[regions[r]] = mat['hbr_opto_post'][:,regions[r],:]/(100-saturation[r]) + 1   
        HbT_wh_preLNAME_all[regions[r]] = mat['hbt_whisker_pre'][:,regions[r],:]/100 + 1   
        HbT_wh_postLNAME_all[regions[r]] = mat['hbt_whisker_post'][:,regions[r],:]/100 + 1       
        HbT_op_preLNAME_all[regions[r]] = mat['hbt_opto_pre'][:,regions[r],:]/100 + 1   
        HbT_op_postLNAME_all[regions[r]] = mat['hbt_opto_post'][:,regions[r],:]/100 + 1   
      
    
    # Extract and shape the arrays for each case
    if area == 'All':   # stack data from all regions into a single 2D array 
        HbO_wh_preLNAME = np.vstack((HbO_wh_preLNAME_all[0], HbO_wh_preLNAME_all[1], HbO_wh_preLNAME_all[2], HbO_wh_preLNAME_all[3]))
        HbO_wh_postLNAME =np.vstack((HbO_wh_postLNAME_all[0], HbO_wh_postLNAME_all[1], HbO_wh_postLNAME_all[2], HbO_wh_postLNAME_all[3]))
        HbO_op_preLNAME = np.vstack((HbO_op_preLNAME_all[0], HbO_op_preLNAME_all[1], HbO_op_preLNAME_all[2], HbO_op_preLNAME_all[3]))
        HbO_op_postLNAME = np.vstack((HbO_op_postLNAME_all[0], HbO_op_postLNAME_all[1], HbO_op_postLNAME_all[2], HbO_op_postLNAME_all[3]))
        HbR_wh_preLNAME = np.vstack((HbR_wh_preLNAME_all[0], HbR_wh_preLNAME_all[1], HbR_wh_preLNAME_all[2], HbR_wh_preLNAME_all[3]))
        HbR_wh_postLNAME = np.vstack((HbR_wh_postLNAME_all[0], HbR_wh_postLNAME_all[1], HbR_wh_postLNAME_all[2], HbR_wh_postLNAME_all[3]))
        HbR_op_preLNAME = np.vstack((HbR_op_preLNAME_all[0], HbR_op_preLNAME_all[1], HbR_op_preLNAME_all[2], HbR_op_preLNAME_all[3]))
        HbR_op_postLNAME = np.vstack((HbR_op_postLNAME_all[0], HbR_op_postLNAME_all[1], HbR_op_postLNAME_all[2], HbR_op_postLNAME_all[3]))
        HbT_wh_preLNAME = np.vstack((HbT_wh_preLNAME_all[0], HbT_wh_preLNAME_all[1], HbT_wh_preLNAME_all[2], HbT_wh_preLNAME_all[3]))
        HbT_wh_postLNAME = np.vstack((HbT_wh_postLNAME_all[0], HbT_wh_postLNAME_all[1], HbT_wh_postLNAME_all[2], HbT_wh_postLNAME_all[3]))
        HbT_op_preLNAME = np.vstack((HbT_op_preLNAME_all[0], HbT_op_preLNAME_all[1], HbT_op_preLNAME_all[2], HbT_op_preLNAME_all[3]))
        HbT_op_postLNAME = np.vstack((HbT_op_postLNAME_all[0], HbT_op_postLNAME_all[1], HbT_op_postLNAME_all[2], HbT_op_postLNAME_all[3]))
    else:               # extract data for specific area 
        HbO_wh_preLNAME = HbO_wh_preLNAME_all[regions[area]]
        HbO_wh_postLNAME = HbO_wh_postLNAME_all[regions[area]]  
        HbO_op_preLNAME = HbO_op_preLNAME_all[regions[area]]
        HbO_op_postLNAME = HbO_op_postLNAME_all[regions[area]] 
        HbR_wh_preLNAME = HbR_wh_preLNAME_all[regions[area]]
        HbR_wh_postLNAME = HbR_wh_postLNAME_all[regions[area]]     
        HbR_op_preLNAME = HbR_op_preLNAME_all[regions[area]] 
        HbR_op_postLNAME = HbR_op_postLNAME_all[regions[area]]  
        HbT_wh_preLNAME = HbT_wh_preLNAME_all[regions[area]] 
        HbT_wh_postLNAME = HbT_wh_postLNAME_all[regions[area]]  
        HbT_op_preLNAME = HbT_op_preLNAME_all[regions[area]] 
        HbT_op_postLNAME = HbT_op_postLNAME_all[regions[area]]
    
    # Find the means and standard deviations
    d.HbO_wh_preLNAME_mean = np.mean(HbO_wh_preLNAME, axis=0)          
    d.HbO_wh_preLNAME_std = np.std(HbO_wh_preLNAME, axis=0)             
    d.HbO_wh_postLNAME_mean = np.mean(HbO_wh_postLNAME, axis=0)          
    d.HbO_wh_postLNAME_std = np.std(HbO_wh_postLNAME, axis=0)   
    d.HbO_op_preLNAME_mean = np.mean(HbO_op_preLNAME, axis=0)          
    d.HbO_op_preLNAME_std = np.std(HbO_op_preLNAME, axis=0)   
    d.HbO_op_postLNAME_mean = np.mean(HbO_op_postLNAME, axis=0)          
    d.HbO_op_postLNAME_std = np.std(HbO_op_postLNAME, axis=0) 
    
    d.HbR_wh_preLNAME_mean = np.mean(HbR_wh_preLNAME, axis=0)          
    d.HbR_wh_preLNAME_std = np.std(HbR_wh_preLNAME, axis=0)             
    d.HbR_wh_postLNAME_mean = np.mean(HbR_wh_postLNAME, axis=0)          
    d.HbR_wh_postLNAME_std = np.std(HbR_wh_postLNAME, axis=0)   
    d.HbR_op_preLNAME_mean = np.mean(HbR_op_preLNAME, axis=0)          
    d.HbR_op_preLNAME_std = np.std(HbR_op_preLNAME, axis=0)   
    d.HbR_op_postLNAME_mean = np.mean(HbR_op_postLNAME, axis=0)          
    d.HbR_op_postLNAME_std = np.std(HbR_op_postLNAME, axis=0)  
    
    d.HbT_wh_preLNAME_mean = np.mean(HbT_wh_preLNAME, axis=0)          
    d.HbT_wh_preLNAME_std = np.std(HbT_wh_preLNAME, axis=0)             
    d.HbT_wh_postLNAME_mean = np.mean(HbT_wh_postLNAME, axis=0)          
    d.HbT_wh_postLNAME_std = np.std(HbT_wh_postLNAME, axis=0)   
    d.HbT_op_preLNAME_mean = np.mean(HbT_op_preLNAME, axis=0)          
    d.HbT_op_preLNAME_std = np.std(HbT_op_preLNAME, axis=0)   
    d.HbT_op_postLNAME_mean = np.mean(HbT_op_postLNAME, axis=0)          
    d.HbT_op_postLNAME_std = np.std(HbT_op_postLNAME, axis=0)  
    
    # Find the standard errors      
    d.HbO_wh_preLNAME_stderr = stats.sem(HbO_wh_preLNAME)                     
    d.HbO_wh_postLNAME_stderr = stats.sem(HbO_wh_postLNAME)           
    d.HbO_op_preLNAME_stderr = stats.sem(HbO_op_preLNAME)           
    d.HbO_op_postLNAME_stderr = stats.sem(HbO_op_postLNAME) 
              
    d.HbR_wh_preLNAME_stderr = stats.sem(HbR_wh_preLNAME)                     
    d.HbR_wh_postLNAME_stderr = stats.sem(HbR_wh_postLNAME)          
    d.HbR_op_preLNAME_stderr = stats.sem(HbR_op_preLNAME)           
    d.HbR_op_postLNAME_stderr = stats.sem(HbR_op_postLNAME)  
      
    d.HbT_wh_preLNAME_stderr = stats.sem(HbT_wh_preLNAME)                    
    d.HbT_wh_postLNAME_stderr = stats.sem(HbT_wh_postLNAME)           
    d.HbT_op_preLNAME_stderr = stats.sem(HbT_op_preLNAME)            
    d.HbT_op_postLNAME_stderr = stats.sem(HbT_op_postLNAME)  

    fig = plt.figure(figsize=(12,10), dpi=80, edgecolor='k') 
    fig.suptitle("LNAME data ({0:s})".format(area))
    
    # Find ylims based on max/min error bar for either graph
    ylim1 = min(min(d.HbR_wh_preLNAME_mean - d.HbR_wh_preLNAME_std), min(d.HbR_wh_postLNAME_mean - d.HbR_wh_postLNAME_std), min(d.HbR_op_preLNAME_mean - d.HbR_op_preLNAME_std), min(d.HbR_op_postLNAME_mean - d.HbR_op_postLNAME_std)) - 0.002
    ylim2 = max(max(d.HbO_wh_preLNAME_mean + d.HbO_wh_preLNAME_std), max(d.HbO_wh_postLNAME_mean + d.HbO_wh_postLNAME_std), max(d.HbO_op_preLNAME_mean + d.HbO_op_preLNAME_std), max(d.HbO_op_postLNAME_mean + d.HbO_op_postLNAME_std)) + 0.002
    
    ax = fig.add_subplot(2,2,1)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Whisker Response - Pre LNAME')  
    plt.errorbar(d.time_LNAME, d.HbO_wh_preLNAME_mean, d.HbO_wh_preLNAME_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time_LNAME, d.HbR_wh_preLNAME_mean, d.HbR_wh_preLNAME_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time_LNAME, d.HbT_wh_preLNAME_mean, d.HbT_wh_preLNAME_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.ylim(ylim1, ylim2)
    plt.legend()
    
    ax = fig.add_subplot(2,2,2)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Whisker Response - Post LNAME')  
    plt.errorbar(d.time_LNAME, d.HbO_wh_postLNAME_mean, d.HbO_wh_postLNAME_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time_LNAME, d.HbR_wh_postLNAME_mean, d.HbR_wh_postLNAME_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time_LNAME, d.HbT_wh_postLNAME_mean, d.HbT_wh_postLNAME_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.ylim(ylim1, ylim2)
    plt.legend()   
    
    ax = fig.add_subplot(2,2,3)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Opto Response - Pre LNAME')  
    plt.errorbar(d.time_LNAME, d.HbO_op_preLNAME_mean, d.HbO_op_preLNAME_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time_LNAME, d.HbR_op_preLNAME_mean, d.HbR_op_preLNAME_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time_LNAME, d.HbT_op_preLNAME_mean, d.HbT_op_preLNAME_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.ylim(ylim1, ylim2)
    plt.legend()  

    ax = fig.add_subplot(2,2,4)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Opto Response - Post LNAME')  
    plt.errorbar(d.time_LNAME, d.HbO_op_postLNAME_mean, d.HbO_op_postLNAME_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time_LNAME, d.HbR_op_postLNAME_mean, d.HbR_op_postLNAME_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time_LNAME, d.HbT_op_postLNAME_mean, d.HbT_op_postLNAME_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.ylim(ylim1, ylim2)
    plt.legend()  
    
    return fig, d

'''Extract the experimental data from the LNAME/7NI data that corresponds to the current model simulation'''
# TODO: change to match new data
# Returns the 3 time variables and 3 hemodynamic variables
# d: class containing the experimental data 
def set_corresponding_LNAME_7NI_experimental_data(d):
    
    time_LNAME = d.time_LNAME 
    
    if p.NeuronType == 'whisker' and p.NOswitch == 'normal':    
        exp_HbO = d.HbO_wh_preLNAME_mean
        exp_HbR = d.HbR_wh_preLNAME_mean
        exp_HbT = d.HbT_wh_preLNAME_mean  
        exp_HbO_std = d.HbO_wh_preLNAME_std
        exp_HbR_std = d.HbR_wh_preLNAME_std 
        exp_HbT_std = d.HbT_wh_preLNAME_std 
    elif p.NeuronType == 'whisker' and p.NOswitch == 'LNAME': 
        exp_HbO = d.HbO_wh_postLNAME_mean
        exp_HbR = d.HbR_wh_postLNAME_mean
        exp_HbT = d.HbT_wh_postLNAME_mean 
        exp_HbO_std = d.HbO_wh_postLNAME_std
        exp_HbR_std = d.HbR_wh_postLNAME_std
        exp_HbT_std = d.HbT_wh_postLNAME_std  
    elif p.NeuronType == 'opto' and p.NOswitch == 'normal':    
        exp_HbO = d.HbO_op_preLNAME_mean
        exp_HbR = d.HbR_op_preLNAME_mean
        exp_HbT = d.HbT_op_preLNAME_mean
        exp_HbO_std = d.HbO_op_preLNAME_std
        exp_HbR_std = d.HbR_op_preLNAME_std
        exp_HbT_std = d.HbT_op_preLNAME_std 
    elif p.NeuronType == 'opto' and p.NOswitch == 'LNAME':
        exp_HbO = d.HbO_op_postLNAME_mean
        exp_HbR = d.HbR_op_postLNAME_mean
        exp_HbT = d.HbT_op_postLNAME_mean  
        exp_HbO_std = d.HbO_op_postLNAME_std
        exp_HbR_std = d.HbR_op_postLNAME_std
        exp_HbT_std = d.HbT_op_postLNAME_std      
#    if p.NeuronType == 'whisker' and p.NOswitch == 'normal':
#        time_HbO = d.time_HbO_wh_preLNAME
#        time_HbR = d.time_HbR_wh_preLNAME
#        time_HbT = d.time_HbT_wh_preLNAME   
#        exp_HbO = d.HbO_wh_preLNAME
#        exp_HbR = d.HbR_wh_preLNAME
#        exp_HbT = d.HbT_wh_preLNAME  
#    elif p.NeuronType == 'whisker' and p.NOswitch == '7NI':
#        time_HbO = d.time_HbO_wh_post7ni
#        time_HbR = d.time_HbR_wh_post7ni
#        time_HbT = d.time_HbT_wh_post7ni   
#        exp_HbO = d.HbO_wh_post7ni
#        exp_HbR = d.HbR_wh_post7ni
#        exp_HbT = d.HbT_wh_post7ni
#    elif p.NeuronType == 'whisker' and p.NOswitch == 'LNAME':
#        time_HbO = d.time_HbO_wh_postLNAME
#        time_HbR = d.time_HbR_wh_postLNAME
#        time_HbT = d.time_HbT_wh_postLNAME   
#        exp_HbO = d.HbO_wh_postLNAME
#        exp_HbR = d.HbR_wh_postLNAME
#        exp_HbT = d.HbT_wh_postLNAME    
#    elif p.NeuronType == 'opto' and p.NOswitch == 'normal':
#        time_HbO = d.time_HbO_op_preLNAME
#        time_HbR = d.time_HbR_op_preLNAME
#        time_HbT = d.time_HbT_op_preLNAME   
#        exp_HbO = d.HbO_op_preLNAME
#        exp_HbR = d.HbR_op_preLNAME
#        exp_HbT = d.HbT_op_preLNAME  
#    elif p.NeuronType == 'opto' and p.NOswitch == '7NI':
#        time_HbO = d.time_HbO_op_post7ni
#        time_HbR = d.time_HbR_op_post7ni
#        time_HbT = d.time_HbT_op_post7ni   
#        exp_HbO = d.HbO_op_post7ni
#        exp_HbR = d.HbR_op_post7ni
#        exp_HbT = d.HbT_op_post7ni
#    elif p.NeuronType == 'opto' and p.NOswitch == 'LNAME':
#        time_HbO = d.time_HbO_op_postLNAME
#        time_HbR = d.time_HbR_op_postLNAME
#        time_HbT = d.time_HbT_op_postLNAME   
#        exp_HbO = d.HbO_op_postLNAME
#        exp_HbR = d.HbR_op_postLNAME
#        exp_HbT = d.HbT_op_postLNAME  
        
    return time_LNAME, exp_HbO, exp_HbR, exp_HbT, exp_HbO_std, exp_HbR_std, exp_HbT_std



















'''Import the NEURAL data from the Berwick experiments using both air and oxygen'''
# Note that the data is averaged over all animals, trials and channels!!
# Returns:
# fig_neural_airOxyComparison_sidebyside - figure of the neural data for air and oxygen experiments in subplots side by side
# fig_neural_airOxyComparison - figure of the neural data for air and oxygen experiments on the same axis
# d - class containing: time_airOxy, neural_air_mean, neural_air_std, neural_oxy_mean, neural_oxy_std
def import_Berwick_Neural_AirOxyData():
    
    mat = io.loadmat('./mat_files/MUA_n13_air_oxygen.mat')
    
    # 13 = number of animals
    # 30 = number of trials for each animal
    # 12 = channels (depth through the cortex - channel 1 is surface of the brain, channel 12 is 1.2mm into the cortex)
    # 250 = data length (time)
       
    # Initialise class to contain the data
    d = Data()
    
    d.time_airOxy = mat['tim_mu'].flatten()               # time vector
    
    '''Find the mean and standard deviation and plot to figure'''

    neural_air_array = mat['wt_2s_air_mua_tot']                  # extract the data from the mat file into a 4D array         
    neural_oxy_array = mat['wt_2s_oxy_mua_tot']

    neural_air_array = np.reshape(neural_air_array, (13*30*12, 250))  # reshape into a 2D array so that the mean and std can be calculated
    neural_oxy_array = np.reshape(neural_oxy_array, (13*30*12, 250))   
    
    d.neural_air_mean = np.mean(neural_air_array, axis=0)           # find the mean profile
    d.neural_air_std = np.std(neural_air_array, axis=0)             # find the standard deviation for all the profiles
    d.neural_oxy_mean = np.mean(neural_oxy_array, axis=0)
    d.neural_oxy_std = np.std(neural_oxy_array, axis=0)
    
    d.neural_air_stderr = stats.sem(neural_air_array)  
    d.neural_oxy_stderr = stats.sem(neural_oxy_array)
    
    # Plot them side by side
    fig_neural_airOxyComparison_sidebyside = plt.figure(dpi=80, edgecolor='k') 
    fig_neural_airOxyComparison_sidebyside.suptitle("Berwick Neural data")
    
    ax = fig_neural_airOxyComparison_sidebyside.add_subplot(1,2,1)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Air')  
    plt.errorbar(d.time_airOxy, d.neural_air_mean, d.neural_air_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b')
    plt.ylim(45,125)
    plt.xlabel('Time [s]') 

    ax = fig_neural_airOxyComparison_sidebyside.add_subplot(1,2,2)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Oxy')  
    plt.errorbar(d.time_airOxy, d.neural_oxy_mean, d.neural_oxy_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r')
    plt.ylim(45,125)
    plt.xlabel('Time [s]')  


    # Plot them both on the same axis
    fig_neural_airOxyComparison = plt.figure(dpi=80, edgecolor='k') 
    fig_neural_airOxyComparison.suptitle("Berwick Neural data")
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) 
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    plt.errorbar(d.time_airOxy, d.neural_air_mean, d.neural_air_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'Air')
    plt.errorbar(d.time_airOxy, d.neural_oxy_mean, d.neural_oxy_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'Oxy')
    plt.xlabel('Time [s]') 
    plt.legend()
    
    return fig_neural_airOxyComparison_sidebyside, fig_neural_airOxyComparison, d




'''Import the data from the new Berwick experiments using both air and oxygen'''
# Returns:
# fig_hemo_airOxyComparison - figure of the averaged hemodynamics for air and oxygen experiments
# d - class containing: time_airOxy, HbOair_mean, HbOair_std, HbRair_mean, HbRair_std, HbTair_mean, HbTair_std,
#                       HbOoxy_mean, HbOoxy_std, HbRoxy_mean, HbRoxy_std, HbToxy_mean, HbToxy_std
# area - choose from: 'All', 'Whisker', 'Artery', 'Vein', 'Parenchyma'. Default = 'All'
def import_Berwick_hemo_AirOxyData(area='All'):
    
    mat = io.loadmat('./mat_files/haems_n13_2s_air_oxygen.mat')
    
    # 13 = number of animals
    # 4 = regions from where the responses are from (1 = whisker region, 2 = Artery, 3 = Vein, 4 = Parenchyma)
    # 200 = data length (time)
    # 30 = the number of trials for each animal
       
    # Initialise class to contain the data
    d = Data()
    
    d.time_airOxy = mat['tim'].flatten() - 5                # time vector
    
    '''Find the mean and standard deviation and plot to figure'''
    
    regions = { 'Whisker': 0, 'Artery': 1, 'Vein': 2, 'Parenchyma': 3} # Possible regions
    
    HbOair_array = mat['wt_hbo_tot_2s_air_Fract']                  # extract the data from the mat file into a 4D array         
    HbRair_array = mat['wt_hbr_tot_2s_air_Fract']
    HbTair_array = mat['wt_hbt_tot_2s_air_Fract']
    HbOoxy_array = mat['wt_hbo_tot_2s_oxy_Fract']
    HbRoxy_array = mat['wt_hbr_tot_2s_oxy_Fract']
    HbToxy_array = mat['wt_hbt_tot_2s_oxy_Fract']
    
    if area == 'All':
        HbOair_array = HbOair_array.transpose(3,1,0,2)          # rearrange data so that it can be reshaped in the next step
        HbOair_array = np.reshape(HbOair_array, (13*4*30, 200)) # reshape into a 2D array so that the mean and std can be calculated
        HbRair_array = HbRair_array.transpose(3,1,0,2)
        HbRair_array = np.reshape(HbRair_array, (13*4*30, 200))   
        HbTair_array = HbTair_array.transpose(3,1,0,2)
        HbTair_array = np.reshape(HbTair_array, (13*4*30, 200))   
        HbOoxy_array = HbOoxy_array.transpose(3,1,0,2)
        HbOoxy_array = np.reshape(HbOoxy_array, (13*4*30, 200))         
        HbRoxy_array = HbRoxy_array.transpose(3,1,0,2)
        HbRoxy_array = np.reshape(HbRoxy_array, (13*4*30, 200))   
        HbToxy_array = HbToxy_array.transpose(3,1,0,2)
        HbToxy_array = np.reshape(HbToxy_array, (13*4*30, 200))   
    else:
        HbOair_array = HbOair_array[:,regions[area],:,:].transpose(2,0,1)   # extract data for specific area then rearrange data so that it can be reshaped in the next step
        HbOair_array = np.reshape(HbOair_array, (13*30, 200))               # reshape into a 2D array so that the mean and std can be calculated
        HbRair_array = HbRair_array[:,regions[area],:,:].transpose(2,0,1)   
        HbRair_array = np.reshape(HbRair_array, (13*30, 200))              
        HbTair_array = HbTair_array[:,regions[area],:,:].transpose(2,0,1)   
        HbTair_array = np.reshape(HbTair_array, (13*30, 200))               
        HbOoxy_array = HbOoxy_array[:,regions[area],:,:].transpose(2,0,1)  
        HbOoxy_array = np.reshape(HbOoxy_array, (13*30, 200))               
        HbRoxy_array = HbRoxy_array[:,regions[area],:,:].transpose(2,0,1)   
        HbRoxy_array = np.reshape(HbRoxy_array, (13*30, 200))              
        HbToxy_array = HbToxy_array[:,regions[area],:,:].transpose(2,0,1)  
        HbToxy_array = np.reshape(HbToxy_array, (13*30, 200))               
    
    d.HbOair_mean = np.mean(HbOair_array, axis=0)           # find the mean Hb profile
    d.HbOair_std = np.std(HbOair_array, axis=0)             # find the standard deviation for all the Hb profiles
    d.HbRair_mean = np.mean(HbRair_array, axis=0)
    d.HbRair_std = np.std(HbRair_array, axis=0)
    d.HbTair_mean = np.mean(HbTair_array, axis=0)
    d.HbTair_std = np.std(HbTair_array, axis=0)
    d.HbOoxy_mean = np.mean(HbOoxy_array, axis=0)
    d.HbOoxy_std = np.std(HbOoxy_array, axis=0)
    d.HbRoxy_mean = np.mean(HbRoxy_array, axis=0)
    d.HbRoxy_std = np.std(HbRoxy_array, axis=0)
    d.HbToxy_mean = np.mean(HbToxy_array, axis=0)
    d.HbToxy_std = np.std(HbToxy_array, axis=0)
    
    # Find the standard error
    d.HbOair_stderr = stats.sem(HbOair_array)             
    d.HbRair_stderr = stats.sem(HbRair_array)
    d.HbTair_stderr = stats.sem(HbTair_array)
    d.HbOoxy_stderr = stats.sem(HbOoxy_array)
    d.HbRoxy_stderr = stats.sem(HbRoxy_array)
    d.HbToxy_stderr = stats.sem(HbToxy_array)
    
    fig_hemo_airOxyComparison = plt.figure(figsize=(12,5), dpi=80, edgecolor='k') 
    fig_hemo_airOxyComparison.suptitle("Berwick data - Haemodynamics (Area: {0:s})".format(area))
    
    ylim1 = 0.96
    ylim2 = 1.07
    
    ax = fig_hemo_airOxyComparison.add_subplot(1,2,1)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Air')  
    plt.errorbar(d.time_airOxy, d.HbOair_mean, d.HbOair_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time_airOxy, d.HbRair_mean, d.HbRair_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time_airOxy, d.HbTair_mean, d.HbTair_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.ylim(ylim1, ylim2)
    plt.legend()
    
    ax = fig_hemo_airOxyComparison.add_subplot(1,2,2)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Oxygen')  
    plt.errorbar(d.time_airOxy, d.HbOoxy_mean, d.HbOoxy_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time_airOxy, d.HbRoxy_mean, d.HbRoxy_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time_airOxy, d.HbToxy_mean, d.HbToxy_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.ylim(ylim1, ylim2)
    plt.legend()
    
    
    return fig_hemo_airOxyComparison, d















'''Import the data from the Berwick experiments for whisker/opto using LNAME, HET'''
# Returns:
# fig_hemo_whiskerOptoComparison - figure of the hemodynamics for whisker and opto experiments
# d - class containing: HbOopto_mean, HbOopto_std, HbRopto_mean, HbRopto_std, HbTopto_mean, HbTopto_std, 
#                       HbOwhisk_mean, HbOwhisk_std, HbRwhisk_mean, HbRwhisk_std, HbTwhisk_mean, HbTwhisk_std, time

# Input parameters:
# dataset - choose from: 'tots_HET_post', 'tots_HET_pre', 'tots_LNAME_HET_post', 'tots_LNAME_HET_pre',
#                        'tots_LNAME_NEURAL_post', 'tots_LNAME_NEURAL_pre', 'tots_LNAME_NON_NEURAL_post', 'tots_LNAME_NON_NEURAL_pre'
#                        'tots_LNAME_post', 'tots_LNAME_pre', 'tots_NODRUG_post', 'tots_NODRUG_pre'
# area - choose from: 'All', 'Whisker', 'Artery', 'Vein', 'Parenchyma'. Default = 'All'
def import_Berwick_HET_LNAME_Data(dataset, area='All'):

    mat = io.loadmat('./mat_files/'+ dataset + '.mat')
    
    # 4 = regions from where the responses are from (1 = whisker region, 2 = Artery, 3 = Vein, 4 = Parenchyma)
    # 200 = data length (time)
       
    # Initialise class to contain the data
    d = Data()
    
    d.time = np.linspace(-5, 20, 200) # Create time vector with 200 timesteps between -5 and 20 sec, stimulation occurring at 0 sec
    
    regions = { 'Whisker': 0, 'Artery': 1, 'Vein': 2, 'Parenchyma': 3} # Possible regions
#    saturation = { 'Whisker': 70, 'Artery': 80, 'Vein': 60, 'Parenchyma': 70} # saturation % for each region, used for scaling the micromolar changes into fractional changes

    
    HbOopto_array = [0]*4 # Initialise the arrays      
    HbRopto_array = [0]*4 
    HbTopto_array = [0]*4 
    HbOwhisk_array = [0]*4 
    HbRwhisk_array = [0]*4 
    HbTwhisk_array = [0]*4 
    
    # Extract and scale data for each region (percentage change from 0 -> percentage change from 1)
    for r in regions:
#        HbOopto_array[regions[r]] = mat['hbo_opto'][:,regions[r],:]/saturation[r] + 1 
#        HbOwhisk_array[regions[r]] = mat['hbo_whisk'][:,regions[r],:]/saturation[r] + 1     
#        HbRopto_array[regions[r]] = mat['hbr_opto'][:,regions[r],:]/(100-saturation[r]) + 1   
#        HbRwhisk_array[regions[r]] = mat['hbr_whisk'][:,regions[r],:]/(100-saturation[r]) + 1       
#        HbTopto_array[regions[r]] = mat['hbt_opto'][:,regions[r],:]/100 + 1   
#        HbTwhisk_array[regions[r]] = mat['hbt_whisk'][:,regions[r],:]/100 + 1       
        
        HbOopto_array[regions[r]] = mat['hbo_opto'][:,regions[r],:]/100 + 1 
        HbOwhisk_array[regions[r]] = mat['hbo_whisk'][:,regions[r],:]/100 + 1     
        HbRopto_array[regions[r]] = mat['hbr_opto'][:,regions[r],:]/100 + 1   
        HbRwhisk_array[regions[r]] = mat['hbr_whisk'][:,regions[r],:]/100 + 1       
        HbTopto_array[regions[r]] = mat['hbt_opto'][:,regions[r],:]/100 + 1   
        HbTwhisk_array[regions[r]] = mat['hbt_whisk'][:,regions[r],:]/100 + 1   
        
        
    # Extract and shape the arrays for each case
    if area == 'All':   # stack data from all regions into a single 2D array 
        HbOopto_array = np.vstack((HbOopto_array[0], HbOopto_array[1], HbOopto_array[2], HbOopto_array[3]))
        HbRopto_array = np.vstack((HbRopto_array[0], HbRopto_array[1], HbRopto_array[2], HbRopto_array[3]))
        HbTopto_array = np.vstack((HbTopto_array[0], HbTopto_array[1], HbTopto_array[2], HbTopto_array[3]))
        HbOwhisk_array = np.vstack((HbOwhisk_array[0], HbOwhisk_array[1], HbOwhisk_array[2], HbOwhisk_array[3]))
        HbRwhisk_array = np.vstack((HbRwhisk_array[0], HbRwhisk_array[1], HbRwhisk_array[2], HbRwhisk_array[3]))
        HbTwhisk_array = np.vstack((HbTwhisk_array[0], HbTwhisk_array[1], HbTwhisk_array[2], HbTwhisk_array[3]))
    else:               # extract data for specific area 
        HbOopto_array = HbOopto_array[regions[area]]
        HbRopto_array = HbRopto_array[regions[area]]
        HbTopto_array = HbTopto_array[regions[area]]
        HbOwhisk_array = HbOwhisk_array[regions[area]]
        HbRwhisk_array = HbRwhisk_array[regions[area]]
        HbTwhisk_array = HbTwhisk_array[regions[area]]

    '''Find the mean and standard deviation'''  
    d.HbOopto_mean = np.mean(HbOopto_array, axis=0)           # find the mean Hb profile
    d.HbOopto_std = np.std(HbOopto_array, axis=0)             # find the standard deviation for all the Hb profiles
    d.HbRopto_mean = np.mean(HbRopto_array, axis=0) 
    d.HbRopto_std = np.std(HbRopto_array, axis=0)
    d.HbTopto_mean = np.mean(HbTopto_array, axis=0) 
    d.HbTopto_std = np.std(HbTopto_array, axis=0)
    d.HbOwhisk_mean = np.mean(HbOwhisk_array, axis=0) 
    d.HbOwhisk_std = np.std(HbOwhisk_array, axis=0)
    d.HbRwhisk_mean = np.mean(HbRwhisk_array, axis=0) 
    d.HbRwhisk_std = np.std(HbRwhisk_array, axis=0)
    d.HbTwhisk_mean = np.mean(HbTwhisk_array, axis=0) 
    d.HbTwhisk_std = np.std(HbTwhisk_array, axis=0)
    
    '''Find the standard error'''
    d.HbOopto_stderr = stats.sem(HbOopto_array)             # find the standard error for all the Hb profiles
    d.HbRopto_stderr = stats.sem(HbRopto_array)
    d.HbTopto_stderr = stats.sem(HbTopto_array)
    d.HbOwhisk_stderr = stats.sem(HbOwhisk_array)
    d.HbRwhisk_stderr = stats.sem(HbRwhisk_array)
    d.HbTwhisk_stderr = stats.sem(HbTwhisk_array)
    
    '''Plot the mean and standard error to figures (Whisker and Opto)'''
    fig_hemo_whiskerOptoComparison = plt.figure(figsize=(12,5), dpi=80, edgecolor='k') 
    fig_hemo_whiskerOptoComparison.suptitle("Berwick data haemodynamics (Dataset: {0:s}, Area: {1:s})".format(dataset, area))
   
    ax = fig_hemo_whiskerOptoComparison.add_subplot(1,2,1)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Whisker')  
    plt.errorbar(d.time, d.HbOwhisk_mean, d.HbOwhisk_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time, d.HbRwhisk_mean, d.HbRwhisk_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time, d.HbTwhisk_mean, d.HbTwhisk_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.legend()
    
    ax = fig_hemo_whiskerOptoComparison.add_subplot(1,2,2)
    ax.get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
    ax.set_title('Opto')  
    plt.errorbar(d.time, d.HbOopto_mean, d.HbOopto_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='r', label = 'HbO')
    plt.errorbar(d.time, d.HbRopto_mean, d.HbRopto_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='b', label = 'HbR')
    plt.errorbar(d.time, d.HbTopto_mean, d.HbTopto_stderr, capsize=3, errorevery=12, elinewidth=1, capthick=1, color='g', label = 'HbT')
    plt.xlabel('Time [s]') 
    plt.legend()

    
    return fig_hemo_whiskerOptoComparison, d
