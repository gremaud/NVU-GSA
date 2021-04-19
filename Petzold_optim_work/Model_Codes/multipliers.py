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

V_IC = np.ones(51)
#data=pd.read_csv('./Samples/V_IC.csv', usecols = [0],names=['name'])
#df = pd.DataFrame(data).astype('float_')
#V_IC = df.values.flatten()
def V_d_function(change_index):
    V_d = np.ones(70+1)
    if change_index >=0:
        V_d[change_index]=0
    
    return V_d



#%% Choose optimisation type
# optimised_once: the top 15 parameters are scaled by a multiplier (obtained by optimising the base model)
# optimised_twice_top15 = the top 15 parameters are scaled by a multiplier (obtained by optimising the 'optimised once' model)
# optimised_twice_all = all parameters are scaled by a mulitplier (obtained by optimising the 'optimised once' model)
# no_optimisation = nominal realisation
def V_c_function(optimisation_type):
    V_c = np.ones(234)
        
    if optimisation_type  == 'optimised_twice_all':
        # Import the multipliers for optimisation from the excel file (optimised twice, all parameters)
        data = pd.read_excel('./September2020ParameterNotes.xlsx')
        df = pd.DataFrame(data, columns = ['Final Optimization Modifier'])
        V_c = df.values.flatten()
    
    elif optimisation_type == 'nominal':
        V_c = np.ones(234)
    
    elif optimisation_type == 'Robin_Pre1':
        locations_matlab=[182, 179, 58, 13, 53, 69, 180, 97, 63, 137, 172,129, 192, 135, 150, 28, 78, 7, 74]
        data=pd.read_csv('./Samples/pre_round1_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]
        
    elif optimisation_type == 'Robin_Pre2':
        
        locations_matlab=[182, 179, 58, 13, 53, 69, 180, 97, 63, 137, 172,129, 192, 135, 150, 28, 78, 7, 74]
        data=pd.read_csv('./Samples/pre_round1_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]    
            
        locations_matlab=[129, 127, 179, 180, 182, 58, 176, 119, 13, 128, 103, 69, 53, 59, 75, 74, 143, 218, 72]
        data=pd.read_csv('./Samples/pre_round2_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]
        
    elif     optimisation_type == 'Robin_Post1':
        
        locations_matlab=[182, 179, 58, 13, 53, 69, 180, 97, 63, 137, 172,129, 192, 135, 150, 28, 78, 7, 74]
        data=pd.read_csv('./Samples/pre_round1_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]    
            
        locations_matlab=[129, 127, 179, 180, 182, 58, 176, 119, 13, 128, 103, 69, 53, 59, 75, 74, 143, 218, 72]
        data=pd.read_csv('./Samples/pre_round2_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]
            
        locations_matlab=[53,58,57,178,179,69,120,63,159]
        data=pd.read_csv('./Samples/post_round1_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]    
    
    elif     optimisation_type == 'Robin_Post2':
        
        locations_matlab=[182, 179, 58, 13, 53, 69, 180, 97, 63, 137, 172,129, 192, 135, 150, 28, 78, 7, 74]
        data=pd.read_csv('./Samples/pre_round1_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]    
            
        locations_matlab=[129, 127, 179, 180, 182, 58, 176, 119, 13, 128, 103, 69, 53, 59, 75, 74, 143, 218, 72]
        data=pd.read_csv('./Samples/pre_round2_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]
            
        locations_matlab=[53,58,57,178,179,69,120,63,159]
        data=pd.read_csv('./Samples/post_round1_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]    
            
        locations_matlab=[180, 179, 63, 182, 58, 53, 212, 178, 69, 61, 67, 143, 70, 144, 102, 137, 56]
        data=pd.read_csv('./Samples/post_round2_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]    
    
    elif optimisation_type == 'optim_in_progress':
        data=pd.read_csv('./Samples/Temp_optim.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data).astype('float_')
        V_c = df.values.flatten()
        
        locations_matlab=[182, 179, 58, 13, 53, 69, 180, 97, 63, 137, 172,129, 192, 135, 150, 28, 78, 7, 74]
        data=pd.read_csv('./Samples/pre_round1_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]    
            
        locations_matlab=[129, 127, 179, 180, 182, 58, 176, 119, 13, 128, 103, 69, 53, 59, 75, 74, 143, 218, 72]
        data=pd.read_csv('./Samples/pre_round2_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]
            
        locations_matlab=[53,58,57,178,179,69,120,63,159]
        data=pd.read_csv('./Samples/post_round1_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]     
        
        
        locations_matlab=[180, 179, 63, 182, 58, 53, 212, 178, 69, 61, 67, 143, 70, 144, 102, 137, 56]
        data=pd.read_csv('./Samples/post_round2_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]
    
    else:
        data=pd.read_csv('./Samples/V_9.csv', usecols = [optimisation_type],names=['name'])
        df = pd.DataFrame(data).astype('float_')
        V_c = df.values.flatten()
        
        locations_matlab=[182, 179, 58, 13, 53, 69, 180, 97, 63, 137, 172,129, 192, 135, 150, 28, 78, 7, 74]
        data=pd.read_csv('./Samples/pre_round1_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]    
            
        locations_matlab=[129, 127, 179, 180, 182, 58, 176, 119, 13, 128, 103, 69, 53, 59, 75, 74, 143, 218, 72]
        data=pd.read_csv('./Samples/pre_round2_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]
            
        locations_matlab=[53,58,57,178,179,69,120,63,159]
        data=pd.read_csv('./Samples/post_round1_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i] 
            
        locations_matlab=[180, 179, 63, 182, 58, 53, 212, 178, 69, 61, 67, 143, 70, 144, 102, 137, 56]
        data=pd.read_csv('./Samples/post_round2_optim_results.csv', usecols = [0],names=['name'])
        df = pd.DataFrame(data)
        temp_values = df.values.flatten()
        for i in range(len(locations_matlab)):
            V_c[locations_matlab[i]-1]=V_c[locations_matlab[i]-1]*temp_values[i]
                         
    return V_c

