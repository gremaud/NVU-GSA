#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 18:14:02 2021

@author: rmmorill
"""

from discrete_sobol_with_holes import discrete_sobol_holes_even_probability as sobol
import numpy as np
import random as rand
import statistics as stats
import matplotlib.pyplot as plt
        
list_of_num_removed=np.arange(10,500,10)

results_absolute=np.zeros(len(list_of_num_removed))
results_realtive_to_mean=np.zeros(len(list_of_num_removed))


for NN in range(len(list_of_num_removed)):
    

    fun_flag=1
    num_vars=10
    num_points_removed=list_of_num_removed[NN]
    
    #num_points_removed=16
    #num_vars=5
    
    if fun_flag==1:
        funf= lambda x: sum(x)
    
    combos=list()
    
    n = num_vars                      
    for i in range(2**n):
        b = bin(i)[2:]
        l = len(b)
        b = str(0) * (n - l) + b
        combos.append(b)
    
    valid_combos=combos.copy()
    for i in range(num_points_removed):
        index=rand.randint(0,len(valid_combos)-1)
        removed=valid_combos[index]
        valid_combos.remove(removed)
    
    values=results_v=np.zeros(len(valid_combos))
    
    for i in range(len(valid_combos)):
        combo=valid_combos[i]
        x=np.zeros(len(combo))
        for j in range(len(combo)):
            if combo[j]=='1':
                x[j]=1
        values[i]=funf(x)
    
    
    
    S , STT =sobol(combos,valid_combos,values)
    
    small=1/num_vars-min(STT)
    large=max(STT)-1/num_vars
    results_absolute[NN]=max([small,large])*num_vars
    
    small=stats.mean(STT)-min(STT)
    large=max(STT)-stats.mean(STT)
    results_realtive_to_mean[NN]=max([small,large])*num_vars
    
    print(num_vars, ' , ', num_points_removed, ' Complete')


fig = plt.figure()
plt.plot(100*list_of_num_removed/(2**num_vars),results_absolute)
plt.xlabel('% points removed') 
plt.ylabel('Relative error (by true ST)') 

fig2 = plt.figure()
plt.plot(100*list_of_num_removed/(2**num_vars),results_realtive_to_mean)
plt.xlabel('% points removed')
plt.ylabel('Relative error (by mean ST)')  