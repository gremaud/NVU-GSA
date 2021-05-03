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
v_a_flag='v'

v_nominal, a_nominal, time =single_eval(values,model,[])
v=[]
a=[]


if v_a_flag =='v':
    for i in range(70+1):
        v_temp, a_temp, time =single_eval(values,model,i)
        v.append(v_temp)
        difference=np.zeros([70+1,2])
    
    for i in range(70+1):
        check=compare_results_v(v[i],v[i],norm_flag)
        if check[0] :
            print('%i failed' % i)
            difference[i]=-1
        else:
            print('%i passed' % i)
            diff=compare_results_v(v_nominal,v[i],norm_flag)
            difference[i,0]=norm(diff[1],norm_flag)
            difference[i,1]=np.argmax(diff[1])

                

if v_a_flag=='a_diff':
    for i in range(70+1):
        v_temp, a_temp, time =single_eval(values,model,i)
        a.append(a_temp)        
    
    
    results=np.zeros([70+1,2])
    for i in range(70+1):
        check=compare_results_a(a[i],a[i],norm_flag)
        if check[0] :
             print('%i failed' % i)
             results[i,:]=-1
        else:
             print('%i passed' % i)
             temp=compare_results_a(a_nominal,a[i],norm_flag)
             results[i,:]=temp[1:]


if v_a_flag == 'a_error':
    if model == 'pre':
        dataset = 'tots_LNAME_pre'
     
    Data = im.import_Berwick_HET_LNAME_Data(dataset, area='Whisker')
    results=np.zeros([70+1,2])
    for i in range(70+1):    
        v_temp, a_temp, time =single_eval(values,model,i)
        check=compare_results_a(a_temp,a_temp,norm_flag)
        if check[0] :
             print('%i failed' % i)
             results[i,:]=-1
        else:
            print('%i passed' % i)
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


            
print('done')


