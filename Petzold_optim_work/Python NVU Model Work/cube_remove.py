import sys
sys.path.append("./Model_Codes/")
from single_evaluation import single_eval
from scipy.interpolate import interp1d
import import_mat_files_no_fig as im

from test_switch_functions import compare_results_v
from test_switch_functions import compare_results_a
    
    
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
v_a_flag='v'


combos=list()
list_of_low_reactions=[5,8,10,15,24,27,30] #Cube 2
#list_of_low_reactions=[4,10,11,12,13,15,17,18,19,22,52,54] # Cube 1


n = len(list_of_low_reactions)                       
for i in range(2**n):
    b = bin(i)[2:]
    l = len(b)
    b = str(0) * (n - l) + b
    combos.append(b)

numpairs=len(combos)
combos_for_ST=combos;

v_nominal, a_nominal, time =single_eval(values,model,[])
if model == 'pre':
    dataset = 'tots_LNAME_pre'
Data = im.import_Berwick_HET_LNAME_Data(dataset, area='Whisker')

if v_a_flag =='v':
    results=np.zeros([numpairs,51])
    QoIs=np.zeros([numpairs])
else:
    results=np.zeros([numpairs,2])

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
    if v_a_flag =='v':
        results[i,:]=compare_results_v(v_nominal,v_temp,inner_norm_flag)
        if min(results[i,:])<0:
            error_num=min(results[i,:])
            QoIs[i]=error_num
            combos_for_ST[i]=error_num
        else:
            QoIs[i]=np.linalg.norm(results[i,:],outer_norm_flag)
    elif v_a_flag=='a_diff':
        results[i,:]=compare_results_a(a_nominal,a_temp,inner_norm_flag)
    elif v_a_flag == 'a_error':
        interpolator=interp1d(time, a_temp.HBO_N)
        HBO_interp=interpolator(Data.time)
        Error_HBO=HBO_interp-Data.HbOwhisk_mean
        Error_HBO=sum(map(lambda x:x*x,Error_HBO))

        interpolator=interp1d(time, a_temp.HBR_N)
        HBR_interp=interpolator(Data.time)
        Error_HBR=HBR_interp-Data.HbRwhisk_mean
        Error_HBR=sum(map(lambda x:x*x,Error_HBR))
        
        results[i,0]=Error_HBO
        results[i,1]=Error_HBR
       


QoIs_for_ST=QoIs;

QoIs_for_ST=[element for element in QoIs_for_ST if element>=0]
combos_for_ST=[element for element in combos_for_ST if not(isinstance(element,float))]
        
ST=[0 for element in list_of_low_reactions]
mu_Total=stats.mean(QoIs_for_ST)
Var_Total=stats.variance(QoIs_for_ST)
for ST_index in range(n):
    ST[ST_index]=0
    temp_combos=list()
    n = len(list_of_low_reactions)-1
    for i in range(2**(n)):
        b = bin(i)[2:]
        l = len(b)
        b = str(0) * (n - l) + b
        temp_combos.append(b)
    
    temp_sum=np.zeros(len(temp_combos))
    for j in range(len(temp_combos)):
        other_digits=temp_combos[j]
        temp1=other_digits[:ST_index] + "1" + other_digits[ST_index:]
        temp0=other_digits[:ST_index] + "0" + other_digits[ST_index:]
        temp_sum[j]=stats.mean([(QoIs_for_ST[combos_for_ST.index(temp1)] if temp1 in combos_for_ST else 0),(QoIs_for_ST[combos_for_ST.index(temp0)] if temp0 in combos_for_ST else 0)])**2


    ST[ST_index]=1-(stats.mean(temp_sum)-mu_Total**2)/Var_Total
    
    
print(list_of_low_reactions)
print(ST)    



