#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 11:59:59 2021

@author: rmmorill
"""



def discrete_sobol_holes_even_probability(combos,valid_combos,values):
    import numpy as np
    import statistics as stats
    p=1/len(valid_combos)
    
    S=np.zeros(len(combos))
    STT=np.zeros(len(combos[0]))
    
    #calculate all the numerators first
    for i in range(len(combos)):
        current_subscript=combos[i]
        S[i]=recurrsive_sum_of_squares_of_sum(current_subscript,combos,valid_combos,values,0)
        S[i]=p * S[i]
        for j in range(len(combos)):
            if j!=i:
                subset_flag=1
                possible_subset=combos[j]
                for index in range(len(possible_subset)):
                    if int(possible_subset[index])>int(current_subscript[index]):
                        subset_flag=0
                        break
                if subset_flag==1:
                    S[i]=S[i]-S[j]
    
    for i in range(len(STT)):
        for j in range(len(S)):
            subset_flag=0
            possible_subset=combos[j]
            if possible_subset[i]=='1':
                STT[i]=STT[i]+S[j]
    
    S=np.delete(S,0)
    Total_var=stats.pvariance(values)
    S=S/Total_var
    STT=STT/Total_var        
    
    return S, STT












def recurrsive_sum_of_squares_of_sum(key,current_set,valid_combos,values,i):
    #if we have done all divisons and are left with a set over which we now
    #need to calculate S^(-1)SUM[F]^2
    if i==(len(key)):
        count=0
        summ=0
        for possible_combo in current_set:
            if possible_combo in valid_combos:
                count=count+1
                summ=summ+values[valid_combos.index(possible_combo)]
        if count==0:
            return 0
        else:
            return (summ**2)/count
    
    #if we might have more divisons necessary
    else:
        #if the i'th variable isn't one to split over, move to next variable
        if key[i]=='0':
            summ=recurrsive_sum_of_squares_of_sum(key,current_set,valid_combos,values,i+1)
        #if the i'th variable is one to split over, create sets and then move 
        #to next varaible in each one
        else:
            set0=list()
            set1=list()
            for combo in current_set:
                if combo[i]=='0':
                    set0.append(combo)
                elif combo[i]=='1':
                    set1.append(combo)
            summ0=recurrsive_sum_of_squares_of_sum(key,set0,valid_combos,values,i+1)
            summ1=recurrsive_sum_of_squares_of_sum(key,set1,valid_combos,values,i+1)
            summ=summ0+summ1
            
        return summ