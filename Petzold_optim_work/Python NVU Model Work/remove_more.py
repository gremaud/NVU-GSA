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

values='Robin_Pre2'
model='pre'

norm_flag=2
v_a_flag='v'

num_removed=2
list_of_all_reactions=range(1,59+1)
combos = list(combinations(list_of_all_reactions, num_removed))


numpairs=len(combos)

v_nominal, a_nominal, time =single_eval(values,model,[])
if model == 'pre':
    dataset = 'tots_LNAME_pre'
Data = im.import_Berwick_HET_LNAME_Data(dataset, area='Whisker')

if v_a_flag =='v':
    results=np.zeros([numpairs,51])
else:
    results=np.zeros([numpairs,2])

for i in range(numpairs):
    print('starting ', i+1, 'of ', numpairs)
    current_removed=list(combos[i])
    v_temp, a_temp, time =single_eval(values,model,current_removed)
    if v_a_flag =='v':
        results[i,:]=compare_results_v(v_nominal,v_temp,norm_flag)
    elif v_a_flag=='a_diff':
        results[i,:]=compare_results_a(a_nominal,a_temp,norm_flag)
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
        