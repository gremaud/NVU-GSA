from single_evaluation import single_eval

def compare_results(v1,v2):
    from numpy.linalg import norm
    diff=np.zeros(51)
    # Neuron
    diff[1-1]=norm(v1.E_t-v2.E_t)/norm(v1.E_t)
    diff[2-1]=norm(v1.I_t-v2.I_t)/norm(v1.I_t)    
    diff[3-1]=norm(v1.K_e-v2.K_e)/norm(v1.K_e)
    diff[4-1]=norm(v1.Na_sa-v2.Na_sa)/norm(v1.Na_sa)
    diff[5-1]=norm(v1.Na_d-v2.Na_d)/norm(v1.Na_d)       
    diff[6-1]=norm(v1.O2-v2.O2)/norm(v1.O2)
    diff[7-1]=norm(v1.CBV-v2.CBV)/norm(v1.CBV)
    diff[8-1]=norm(v1.HbR-v2.HbR)/norm(v1.HbR)
    diff[9-1]=norm(v1.Ca_n-v2.Ca_n)/norm(v1.Ca_n)
    diff[10-1]=norm(v1.nNOS-v2.nNOS)/norm(v1.nNOS) 
    diff[11-1]=norm(v1.NO_n-v2.NO_n)/norm(v1.NO_n)
    
    # Astrocyte
    diff[12-1]=norm(v1.v_k-v2.v_k)/norm(v1.v_k)
    diff[13-1]=norm(v1.K_p-v2.K_p)/norm(v1.K_p)
    diff[14-1]=norm(v1.Ca_p-v2.Ca_p)/norm(v1.Ca_p)
    diff[15-1]=norm(v1.Na_k-v2.Na_k)/norm(v1.Na_k)
    diff[16-1]=norm(v1.K_k-v2.K_k)/norm(v1.K_k)
    diff[17-1]=norm(v1.Cl_k-v2.Cl_k)/norm(v1.Cl_k)
    diff[18-1]=norm(v1.HCO3_k-v2.HCO3_k)/norm(v1.HCO3_k)
    diff[19-1]=norm(v1.Na_s-v2.Na_s)/norm(v1.Na_s)
    diff[20-1]=norm(v1.K_s-v2.K_s)/norm(v1.K_s)
    diff[21-1]=norm(v1.HCO3_s-v2.HCO3_s)/norm(v1.HCO3_s)
    diff[22-1]=norm(v1.w_k-v2.w_k)/norm(v1.w_k)
    diff[23-1]=norm(v1.I_k-v2.I_k)/norm(v1.I_k)
    diff[24-1]=norm(v1.Ca_k-v2.Ca_k)/norm(v1.Ca_k)
    diff[25-1]=norm(v1.h_k-v2.h_k)/norm(v1.h_k)
    diff[26-1]=norm(v1.s_k-v2.s_k)/norm(v1.s_k)
    diff[27-1]=norm(v1.m_k-v2.m_k)/norm(v1.m_k)
    diff[28-1]=norm(v1.eet_k-v2.eet_k)/norm(v1.eet_k)
    diff[29-1]=norm(v1.NO_k-v2.NO_k)/norm(v1.NO_k)
    diff[30-1]=norm(v1.AA_k-v2.AA_k)/norm(v1.AA_k)
    
    # SMC
    diff[31-1]=norm(v1.Ca_i-v2.Ca_i)/norm(v1.Ca_i)
    diff[32-1]=norm(v1.s_i-v2.s_i)/norm(v1.s_i)
    diff[33-1]=norm(v1.v_i-v2.v_i)/norm(v1.v_i)
    diff[34-1]=norm(v1.w_i-v2.w_i)/norm(v1.w_i)
    diff[35-1]=norm(v1.I_i-v2.I_i)/norm(v1.I_i)
    diff[36-1]=norm(v1.NO_i-v2.NO_i)/norm(v1.NO_i)
    diff[37-1]=norm(v1.E_b-v2.E_b)/norm(v1.E_b)
    diff[38-1]=norm(v1.E_6c-v2.E_6c)/norm(v1.E_6c)
    diff[39-1]=norm(v1.cGMP_i-v2.cGMP_i)/norm(v1.cGMP_i)
    diff[40-1]=norm(v1.H_i-v2.H_i)/norm(v1.H_i)
    diff[41-1]=norm(v1.AA_i-v2.AA_i)/norm(v1.AA_i)
    
    # EC
    diff[42-1]=norm(v1.Ca_j-v2.Ca_j)/norm(v1.Ca_j)
    diff[43-1]=norm(v1.s_j-v2.s_j)/norm(v1.s_j)
    diff[44-1]=norm(v1.v_j-v2.v_j)/norm(v1.v_j)
    diff[45-1]=norm(v1.I_j-v2.I_j)/norm(v1.I_j)
    diff[46-1]=norm(v1.eNOS-v2.eNOS)/norm(v1.eNOS)
    diff[47-1]=norm(v1.NO_j-v2.NO_j)/norm(v1.NO_j)
    
    # Wall Mechanics
    diff[48-1]=norm(v1.Mp-v2.Mp)/norm(v1.Mp)
    diff[49-1]=norm(v1.AMp-v2.AMp)/norm(v1.AMp)
    diff[50-1]=norm(v1.AM-v2.AM)/norm(v1.AM)
    diff[51-1]=norm(v1.R-v2.R)/norm(v1.R)
    
    NaNcheck=math.isnan(sum(diff))
    mean_diff=100*sum(diff)/51
    max_diff=100*max(diff)
    
    return NaNcheck, mean_diff, max_diff, diff



import math
import numpy as np

values='nominal'
model='pre'

v_nominal=single_eval(values,model,-1)
v=[]
for i in range(70+1):
    v_temp=single_eval(values,model,i)
    v.append(v_temp)

difference=np.zeros([70+1,3])
for i in range(70+1):
    check=compare_results(v[i],v[i])
    if check[0] :
        print('%i failed' % i)
        difference[i,1]=-1
        difference[i,2]=-1
    else:
        print('%i passed' % i)
        difference[i,0], difference[i,1], difference[i,2], diff=compare_results(v_nominal,v[i])
        #print('Mean Difference: %f' % difference[i,1])
        #print('Max Difference: %f' % difference[i,2])
        
print('done')


