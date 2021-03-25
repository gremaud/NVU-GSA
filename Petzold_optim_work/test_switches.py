from single_evaluation import single_eval

def compare_results_v(v1,v2,norm_flag):
    from numpy.linalg import norm
    diff=np.zeros(51)
    # Neuron
    diff[1-1]=norm(v1.E_t-v2.E_t,norm_flag)/norm(v1.E_t,norm_flag)
    diff[2-1]=norm(v1.I_t-v2.I_t,norm_flag)/norm(v1.I_t,norm_flag)    
    diff[3-1]=norm(v1.K_e-v2.K_e,norm_flag)/norm(v1.K_e,norm_flag)
    diff[4-1]=norm(v1.Na_sa-v2.Na_sa,norm_flag)/norm(v1.Na_sa,norm_flag)
    diff[5-1]=norm(v1.Na_d-v2.Na_d,norm_flag)/norm(v1.Na_d,norm_flag)       
    diff[6-1]=norm(v1.O2-v2.O2,norm_flag)/norm(v1.O2,norm_flag)
    diff[7-1]=norm(v1.CBV-v2.CBV,norm_flag)/norm(v1.CBV,norm_flag)
    diff[8-1]=norm(v1.HbR-v2.HbR,norm_flag)/norm(v1.HbR,norm_flag)
    diff[9-1]=norm(v1.Ca_n-v2.Ca_n,norm_flag)/norm(v1.Ca_n,norm_flag)
    diff[10-1]=norm(v1.nNOS-v2.nNOS,norm_flag)/norm(v1.nNOS,norm_flag) 
    diff[11-1]=norm(v1.NO_n-v2.NO_n,norm_flag)/norm(v1.NO_n,norm_flag)
    
    # Astrocyte
    diff[12-1]=norm(v1.v_k-v2.v_k,norm_flag)/norm(v1.v_k,norm_flag)
    diff[13-1]=norm(v1.K_p-v2.K_p,norm_flag)/norm(v1.K_p,norm_flag)
    diff[14-1]=norm(v1.Ca_p-v2.Ca_p,norm_flag)/norm(v1.Ca_p,norm_flag)
    diff[15-1]=norm(v1.Na_k-v2.Na_k,norm_flag)/norm(v1.Na_k,norm_flag)
    diff[16-1]=norm(v1.K_k-v2.K_k,norm_flag)/norm(v1.K_k,norm_flag)
    diff[17-1]=norm(v1.Cl_k-v2.Cl_k,norm_flag)/norm(v1.Cl_k,norm_flag)
    diff[18-1]=norm(v1.HCO3_k-v2.HCO3_k,norm_flag)/norm(v1.HCO3_k,norm_flag)
    diff[19-1]=norm(v1.Na_s-v2.Na_s,norm_flag)/norm(v1.Na_s,norm_flag)
    diff[20-1]=norm(v1.K_s-v2.K_s,norm_flag)/norm(v1.K_s,norm_flag)
    diff[21-1]=norm(v1.HCO3_s-v2.HCO3_s,norm_flag)/norm(v1.HCO3_s,norm_flag)
    diff[22-1]=norm(v1.w_k-v2.w_k,norm_flag)/norm(v1.w_k,norm_flag)
    diff[23-1]=norm(v1.I_k-v2.I_k,norm_flag)/norm(v1.I_k,norm_flag)
    diff[24-1]=norm(v1.Ca_k-v2.Ca_k,norm_flag)/norm(v1.Ca_k,norm_flag)
    diff[25-1]=norm(v1.h_k-v2.h_k,norm_flag)/norm(v1.h_k,norm_flag)
    diff[26-1]=norm(v1.s_k-v2.s_k,norm_flag)/norm(v1.s_k,norm_flag)
    diff[27-1]=norm(v1.m_k-v2.m_k,norm_flag)/norm(v1.m_k,norm_flag)
    diff[28-1]=norm(v1.eet_k-v2.eet_k,norm_flag)/norm(v1.eet_k,norm_flag)
    diff[29-1]=norm(v1.NO_k-v2.NO_k,norm_flag)/norm(v1.NO_k,norm_flag)
    diff[30-1]=norm(v1.AA_k-v2.AA_k,norm_flag)/norm(v1.AA_k,norm_flag)
    
    # SMC
    diff[31-1]=norm(v1.Ca_i-v2.Ca_i,norm_flag)/norm(v1.Ca_i,norm_flag)
    diff[32-1]=norm(v1.s_i-v2.s_i,norm_flag)/norm(v1.s_i,norm_flag)
    diff[33-1]=norm(v1.v_i-v2.v_i,norm_flag)/norm(v1.v_i,norm_flag)
    diff[34-1]=norm(v1.w_i-v2.w_i,norm_flag)/norm(v1.w_i,norm_flag)
    diff[35-1]=norm(v1.I_i-v2.I_i,norm_flag)/norm(v1.I_i,norm_flag)
    diff[36-1]=norm(v1.NO_i-v2.NO_i,norm_flag)/norm(v1.NO_i,norm_flag)
    diff[37-1]=norm(v1.E_b-v2.E_b,norm_flag)/norm(v1.E_b,norm_flag)
    diff[38-1]=norm(v1.E_6c-v2.E_6c,norm_flag)/norm(v1.E_6c,norm_flag)
    diff[39-1]=norm(v1.cGMP_i-v2.cGMP_i,norm_flag)/norm(v1.cGMP_i,norm_flag)
    diff[40-1]=norm(v1.H_i-v2.H_i,norm_flag)/norm(v1.H_i,norm_flag)
    diff[41-1]=norm(v1.AA_i-v2.AA_i,norm_flag)/norm(v1.AA_i,norm_flag)
    
    # EC
    diff[42-1]=norm(v1.Ca_j-v2.Ca_j,norm_flag)/norm(v1.Ca_j,norm_flag)
    diff[43-1]=norm(v1.s_j-v2.s_j,norm_flag)/norm(v1.s_j,norm_flag)
    diff[44-1]=norm(v1.v_j-v2.v_j,norm_flag)/norm(v1.v_j,norm_flag)
    diff[45-1]=norm(v1.I_j-v2.I_j,norm_flag)/norm(v1.I_j,norm_flag)
    diff[46-1]=norm(v1.eNOS-v2.eNOS,norm_flag)/norm(v1.eNOS,norm_flag)
    diff[47-1]=norm(v1.NO_j-v2.NO_j,norm_flag)/norm(v1.NO_j,norm_flag)
    
    # Wall Mechanics
    diff[48-1]=norm(v1.Mp-v2.Mp,norm_flag)/norm(v1.Mp,norm_flag)
    diff[49-1]=norm(v1.AMp-v2.AMp,norm_flag)/norm(v1.AMp,norm_flag)
    diff[50-1]=norm(v1.AM-v2.AM,norm_flag)/norm(v1.AM,norm_flag)
    diff[51-1]=norm(v1.R-v2.R,norm_flag)/norm(v1.R,norm_flag)
    
    NaNcheck=math.isnan(sum(diff))
    
    return NaNcheck, diff

def compare_results_a(a1,a2,norm_flag):
    from numpy.linalg import norm
    import import_mat_files as im
    
    dataset = 'tots_NODRUG_pre'
    fig_hemo_whiskerOptoComparison, LNAME_HET_Data = im.import_Berwick_HET_LNAME_Data(dataset, area='Whisker')

    
    diff_HBO=norm(a1.HBO_N-a2.HBO_N,norm_flag)/norm(a1.HBO_N,norm_flag)
    diff_HBR=norm(a1.HBR_N-a2.HBR_N,norm_flag)/norm(a1.HBR_N,norm_flag)
    
    NaNcheck=math.isnan(diff_HBO + diff_HBR)
    
    
    return [NaNcheck,diff_HBO, diff_HBR]
    
    


import math
import numpy as np
from numpy.linalg import norm

values='Robin_Pre2'
model='pre'

norm_flag=np.inf

v_nominal, a_nominal =single_eval(values,model,-1)
v=[]
a=[]



# for i in range(70+1):
#     v_temp, a_temp =single_eval(values,model,i)
#     v.append(v_temp)
# if norm_flag==np.inf:
#     difference=np.zeros([70+1,2])
# else:
#     difference=np.zeros([70+1])

# for i in range(70+1):
#     check=compare_results_v(v[i],v[i],norm_flag)
#     if check[0] :
#         print('%i failed' % i)
#         difference[i]=-1
#     else:
#         print('%i passed' % i)
#         diff=compare_results_v(v_nominal,v[i],norm_flag)
#         if norm_flag==np.inf:
#             difference[i,0]=norm(diff[1],norm_flag)
#             difference[i,1]=np.argmax(diff[1])
#         else:
#             difference[i]=norm(diff[1],norm_flag)

for i in range(70+1):
    v_temp, a_temp =single_eval(values,model,i)
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

            
print('done')


