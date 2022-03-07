import sys
sys.path.append("./Model_Codes/")
from single_evaluation import single_eval
from scipy.interpolate import interp1d
import import_mat_files_no_fig as im

from test_switch_functions import compare_results_v
from test_switch_functions import compare_results_a
from test_switch_functions import compare_results_sobolev
from discrete_sobol_with_holes import discrete_sobol_holes_even_probability as sobol
    
import math
import numpy as np
from numpy.linalg import norm
from itertools import combinations
import statistics as stats

values='Robin_Pre2'
model='pre'

inner_norm_flag=2
#outer_norm_flag=1
outer_norm_flag=np.inf


combos=list()
#list_of_low_reactions=[5,8,10,15,24,27,30] #Cube 2
#list_of_low_reactions=[4,10,11,12,13,15,17,18,19,22,52,54] # Cube 1
#list_of_low_reactions=[59,58,48,42,4,11]

list_of_low_reactions=[5,8,10]




n = len(list_of_low_reactions)                       
for i in range(2**n):
    b = bin(i)[2:]
    l = len(b)
    b = str(0) * (n - l) + b
    combos.append(b)

numpairs=len(combos)
combos_for_ST=combos.copy();

v_nominal, a_nominal, time =single_eval(values,model,[])
if model == 'pre':
    dataset = 'tots_LNAME_pre'
Data = im.import_Berwick_HET_LNAME_Data(dataset, area='Whisker')


results_v=np.zeros([numpairs,51])
QoIs=np.zeros([numpairs,3])
sobolev_norm_R=np.zeros([numpairs,1])
sobolev_norm_HBR=np.zeros([numpairs,1])
sobolev_norm_Ca_i=np.zeros([numpairs,1])

list_removed=list()
for i in range(numpairs):
    print('starting ', i+1, 'of ', numpairs)
    
    current_removed_index=list(combos[i])
    current_removed=list()
    for j in range(n):
        if current_removed_index[j]=='1':
            current_removed.append(list_of_low_reactions[j])
    list_removed.append(current_removed)
            
    v_temp, a_temp, time =single_eval(values,model,current_removed)

    
        

    interpolator=interp1d(time, a_temp.HBO_N)
    HBO_interp=interpolator(Data.time)
    Error_HBO=HBO_interp-Data.HbOwhisk_mean
    Error_HBO=sum(map(lambda x:x*x,Error_HBO))

    interpolator=interp1d(time, a_temp.HBR_N)
    HBR_interp=interpolator(Data.time)
    Error_HBR=HBR_interp-Data.HbRwhisk_mean
    Error_HBR=sum(map(lambda x:x*x,Error_HBR))
    
    QoIs[i,1]=Error_HBO
    QoIs[i,2]=Error_HBR
    
    results_v[i,:]=compare_results_sobolev(v_nominal,v_temp)
    if min(results_v[i,:])<0:
        error_num=min(results_v[i,:])
        QoIs[i,0]=error_num
        combos_for_ST[i]=error_num
        
        QoIs[i,1]=error_num
        QoIs[i,2]=error_num
    else:
        sobolev_norm_HBR[i]=results_v[i,7]
        sobolev_norm_Ca_i[i]=results_v[i,30]
        sobolev_norm_R[i]=results_v[i,50]
        QoIs[i,0]=sobolev_norm_R[i]+sobolev_norm_HBR[i]+sobolev_norm_Ca_i[i]
       


QoIs_for_ST=QoIs;


QoIs_for_ST0=[element for element in QoIs_for_ST[:,0] if element>=0]
QoIs_for_ST1=[element for element in QoIs_for_ST[:,1] if element>=0]
QoIs_for_ST2=[element for element in QoIs_for_ST[:,2] if element>=0]
combos_for_ST=[element for element in combos_for_ST if not(isinstance(element,float))]

ST=np.zeros([len(list_of_low_reactions),3])

tempS, tempST = sobol(combos,combos_for_ST,QoIs_for_ST0)
ST[:,0]=tempST
tempS, tempST = sobol(combos,combos_for_ST,QoIs_for_ST1)
ST[:,1]=tempST
tempS, tempST = sobol(combos,combos_for_ST,QoIs_for_ST2)
ST[:,2]=tempST
    
print(list_of_low_reactions)
print('Out of ' ,len(combos), ' combinations, ', len(combos_for_ST), ' did not have errors;', len(combos_for_ST)/len(combos)*100,'%')
print(ST[:,0])  
print(ST[:,1])
print(ST[:,2])      



