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

values='Robin_Pre2'
model='pre'

norm_flag=2
v_a_flag='a_error'

#list_of_low_reactions=[34,35,68,29,4,32,63,11,26,31,13,48,8,27,18,28,49,30,65]
list_of_low_reactions=range(1,70+1)
numpairs=int((len(list_of_low_reactions)*len(list_of_low_reactions)-1)/2)

v_nominal, a_nominal, time =single_eval(values,model,[])
if model == 'pre':
    dataset = 'tots_LNAME_pre'
Data = im.import_Berwick_HET_LNAME_Data(dataset, area='Whisker')

results=np.zeros([numpairs,2])
pair_list=[]
index=0
for i in range(len(list_of_low_reactions)):
    for j in range(i+1,len(list_of_low_reactions)):
        pair_list.append([list_of_low_reactions[i],list_of_low_reactions[j]])
        v_temp, a_temp, time =single_eval(values,model,[i,j])
        check_v=compare_results_v(v_temp,v_temp,norm_flag)
        check_a=compare_results_a(a_temp,a_temp,norm_flag)
        if v_a_flag =='v' and check_v[0]:
            results[index,0]=-1
            results[index,1]=-1
            print(pair_list[-1], ' failed')
        elif v_a_flag =='a_diff' and check_a[0]:
            results[index,0]=-1
            results[index,1]=-1
            print(pair_list[-1], ' failed')
        elif v_a_flag =='a_error' and check_a[0]:
            results[index,0]=-1
            results[index,1]=-1
            print(pair_list[-1], ' failed')
        else:
            if v_a_flag =='v':
                diff=compare_results_v(v_nominal,v_temp,norm_flag)
                results[index,0]=norm(diff[1],norm_flag)
                results[index,1]=np.argmax(diff[1])
            elif v_a_flag=='a_diff':
                temp=compare_results_a(a_nominal,a_temp,norm_flag)
                results[index,:]=temp[1:]
            elif v_a_flag == 'a_error':
                interpolator=interp1d(time, a_temp.HBO_N)
                HBO_interp=interpolator(Data.time)
                Error_HBO=HBO_interp-Data.HbOwhisk_mean
                Error_HBO=sum(map(lambda x:x*x,Error_HBO))
    
                interpolator=interp1d(time, a_temp.HBR_N)
                HBR_interp=interpolator(Data.time)
                Error_HBR=HBR_interp-Data.HbRwhisk_mean
                Error_HBR=sum(map(lambda x:x*x,Error_HBR))
                
                results[index,0]=Error_HBO
                results[index,1]=Error_HBR
            print(pair_list[-1], ' passed')
        index+=1

# v_0, a_0, time = single_eval('Robin_Pre2','pre',0)
# check=compare_results_v(v_0,v_0,norm_flag)
# diff=compare_results_v(v_nominal,v_0,norm_flag)


# v_3, a_3, time = single_eval('Robin_Pre2','pre',3)
# v_4, a_4, time = single_eval('Robin_Pre2','pre',4)

# v_03, a_03, time = single_eval('Robin_Pre2','pre',[0,3])
# v_04, a_04, time = single_eval('Robin_Pre2','pre',[0,4])

# check0=compare_results_v(v_nominal,v_0,norm_flag)

# check1=compare_results_v(v_3,v_4,norm_flag)
# check2=compare_results_a(a_nominal,a_3,norm_flag)

# check3=compare_results_v(v_3,v_03,norm_flag)
# check3b=compare_results_a(a_3,a_03,norm_flag)
# check4=compare_results_v(v_4,v_04,norm_flag)
# check4b=compare_results_a(a_4,a_04,norm_flag)