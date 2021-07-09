def compare_results_v(v1,v2,norm_flag):
    import math
    import numpy as np
    from numpy.linalg import norm
    t = np.arange(start=0, stop=550e3+0.001e3, step=0.001e3)
    time=  t/1e3-500e3/1e3
    xlim1=np.where(time>=-5.0)
    xlim1=xlim1[0][0]
    xlim2=np.where(time>=20.001)
    xlim2=xlim2[0][0]
    del t
    diff=np.zeros(51)
    # Neuron
    diff[1-1]=norm(v1.E_t[xlim1:xlim2]-v2.E_t[xlim1:xlim2],norm_flag)/norm(v1.E_t[xlim1:xlim2],norm_flag)
    diff[2-1]=norm(v1.I_t[xlim1:xlim2]-v2.I_t[xlim1:xlim2],norm_flag)/norm(v1.I_t[xlim1:xlim2],norm_flag)    
    diff[3-1]=norm(v1.K_e[xlim1:xlim2]-v2.K_e[xlim1:xlim2],norm_flag)/norm(v1.K_e[xlim1:xlim2],norm_flag)
    diff[4-1]=norm(v1.Na_sa[xlim1:xlim2]-v2.Na_sa[xlim1:xlim2],norm_flag)/norm(v1.Na_sa[xlim1:xlim2],norm_flag)
    diff[5-1]=norm(v1.Na_d[xlim1:xlim2]-v2.Na_d[xlim1:xlim2],norm_flag)/norm(v1.Na_d[xlim1:xlim2],norm_flag)       
    diff[6-1]=norm(v1.O2[xlim1:xlim2]-v2.O2[xlim1:xlim2],norm_flag)/norm(v1.O2[xlim1:xlim2],norm_flag)
    diff[7-1]=norm(v1.CBV[xlim1:xlim2]-v2.CBV[xlim1:xlim2],norm_flag)/norm(v1.CBV[xlim1:xlim2],norm_flag)
    diff[8-1]=norm(v1.HbR[xlim1:xlim2]-v2.HbR[xlim1:xlim2],norm_flag)/norm(v1.HbR[xlim1:xlim2],norm_flag)
    diff[9-1]=norm(v1.Ca_n[xlim1:xlim2]-v2.Ca_n[xlim1:xlim2],norm_flag)/norm(v1.Ca_n[xlim1:xlim2],norm_flag)
    diff[10-1]=norm(v1.nNOS[xlim1:xlim2]-v2.nNOS[xlim1:xlim2],norm_flag)/norm(v1.nNOS[xlim1:xlim2],norm_flag) 
    diff[11-1]=norm(v1.NO_n[xlim1:xlim2]-v2.NO_n[xlim1:xlim2],norm_flag)/norm(v1.NO_n[xlim1:xlim2],norm_flag)
    
    # Astrocyte
    diff[12-1]=norm(v1.v_k[xlim1:xlim2]-v2.v_k[xlim1:xlim2],norm_flag)/norm(v1.v_k[xlim1:xlim2],norm_flag)
    diff[13-1]=norm(v1.K_p[xlim1:xlim2]-v2.K_p[xlim1:xlim2],norm_flag)/norm(v1.K_p[xlim1:xlim2],norm_flag)
    diff[14-1]=norm(v1.Ca_p[xlim1:xlim2]-v2.Ca_p[xlim1:xlim2],norm_flag)/norm(v1.Ca_p[xlim1:xlim2],norm_flag)
    diff[15-1]=norm(v1.Na_k[xlim1:xlim2]-v2.Na_k[xlim1:xlim2],norm_flag)/norm(v1.Na_k[xlim1:xlim2],norm_flag)
    diff[16-1]=norm(v1.K_k[xlim1:xlim2]-v2.K_k[xlim1:xlim2],norm_flag)/norm(v1.K_k[xlim1:xlim2],norm_flag)
    diff[17-1]=norm(v1.Cl_k[xlim1:xlim2]-v2.Cl_k[xlim1:xlim2],norm_flag)/norm(v1.Cl_k[xlim1:xlim2],norm_flag)
    diff[18-1]=norm(v1.HCO3_k[xlim1:xlim2]-v2.HCO3_k[xlim1:xlim2],norm_flag)/norm(v1.HCO3_k[xlim1:xlim2],norm_flag)
    diff[19-1]=norm(v1.Na_s[xlim1:xlim2]-v2.Na_s[xlim1:xlim2],norm_flag)/norm(v1.Na_s[xlim1:xlim2],norm_flag)
    diff[20-1]=norm(v1.K_s[xlim1:xlim2]-v2.K_s[xlim1:xlim2],norm_flag)/norm(v1.K_s[xlim1:xlim2],norm_flag)
    diff[21-1]=norm(v1.HCO3_s[xlim1:xlim2]-v2.HCO3_s[xlim1:xlim2],norm_flag)/norm(v1.HCO3_s[xlim1:xlim2],norm_flag)
    diff[22-1]=norm(v1.w_k[xlim1:xlim2]-v2.w_k[xlim1:xlim2],norm_flag)/norm(v1.w_k[xlim1:xlim2],norm_flag)
    diff[23-1]=norm(v1.I_k[xlim1:xlim2]-v2.I_k[xlim1:xlim2],norm_flag)/norm(v1.I_k[xlim1:xlim2],norm_flag)
    diff[24-1]=norm(v1.Ca_k[xlim1:xlim2]-v2.Ca_k[xlim1:xlim2],norm_flag)/norm(v1.Ca_k[xlim1:xlim2],norm_flag)
    diff[25-1]=norm(v1.h_k[xlim1:xlim2]-v2.h_k[xlim1:xlim2],norm_flag)/norm(v1.h_k[xlim1:xlim2],norm_flag)
    diff[26-1]=norm(v1.s_k[xlim1:xlim2]-v2.s_k[xlim1:xlim2],norm_flag)/norm(v1.s_k[xlim1:xlim2],norm_flag)
    diff[27-1]=norm(v1.m_k[xlim1:xlim2]-v2.m_k[xlim1:xlim2],norm_flag)/norm(v1.m_k[xlim1:xlim2],norm_flag)
    diff[28-1]=norm(v1.eet_k[xlim1:xlim2]-v2.eet_k[xlim1:xlim2],norm_flag)/norm(v1.eet_k[xlim1:xlim2],norm_flag)
    diff[29-1]=norm(v1.NO_k[xlim1:xlim2]-v2.NO_k[xlim1:xlim2],norm_flag)/norm(v1.NO_k[xlim1:xlim2],norm_flag)
    diff[30-1]=norm(v1.AA_k[xlim1:xlim2]-v2.AA_k[xlim1:xlim2],norm_flag)/norm(v1.AA_k[xlim1:xlim2],norm_flag)
    
    # SMC
    diff[31-1]=norm(v1.Ca_i[xlim1:xlim2]-v2.Ca_i[xlim1:xlim2],norm_flag)/norm(v1.Ca_i[xlim1:xlim2],norm_flag)
    diff[32-1]=norm(v1.s_i[xlim1:xlim2]-v2.s_i[xlim1:xlim2],norm_flag)/norm(v1.s_i[xlim1:xlim2],norm_flag)
    diff[33-1]=norm(v1.v_i[xlim1:xlim2]-v2.v_i[xlim1:xlim2],norm_flag)/norm(v1.v_i[xlim1:xlim2],norm_flag)
    diff[34-1]=norm(v1.w_i[xlim1:xlim2]-v2.w_i[xlim1:xlim2],norm_flag)/norm(v1.w_i[xlim1:xlim2],norm_flag)
    diff[35-1]=norm(v1.I_i[xlim1:xlim2]-v2.I_i[xlim1:xlim2],norm_flag)/norm(v1.I_i[xlim1:xlim2],norm_flag)
    diff[36-1]=norm(v1.NO_i[xlim1:xlim2]-v2.NO_i[xlim1:xlim2],norm_flag)/norm(v1.NO_i[xlim1:xlim2],norm_flag)
    diff[37-1]=norm(v1.E_b[xlim1:xlim2]-v2.E_b[xlim1:xlim2],norm_flag)/norm(v1.E_b[xlim1:xlim2],norm_flag)
    diff[38-1]=norm(v1.E_6c[xlim1:xlim2]-v2.E_6c[xlim1:xlim2],norm_flag)/norm(v1.E_6c[xlim1:xlim2],norm_flag)
    diff[39-1]=norm(v1.cGMP_i[xlim1:xlim2]-v2.cGMP_i[xlim1:xlim2],norm_flag)/norm(v1.cGMP_i[xlim1:xlim2],norm_flag)
    diff[40-1]=norm(v1.H_i[xlim1:xlim2]-v2.H_i[xlim1:xlim2],norm_flag)/norm(v1.H_i[xlim1:xlim2],norm_flag)
    diff[41-1]=norm(v1.AA_i[xlim1:xlim2]-v2.AA_i[xlim1:xlim2],norm_flag)/norm(v1.AA_i[xlim1:xlim2],norm_flag)
    
    # EC
    diff[42-1]=norm(v1.Ca_j[xlim1:xlim2]-v2.Ca_j[xlim1:xlim2],norm_flag)/norm(v1.Ca_j[xlim1:xlim2],norm_flag)
    diff[43-1]=norm(v1.s_j[xlim1:xlim2]-v2.s_j[xlim1:xlim2],norm_flag)/norm(v1.s_j[xlim1:xlim2],norm_flag)
    diff[44-1]=norm(v1.v_j[xlim1:xlim2]-v2.v_j[xlim1:xlim2],norm_flag)/norm(v1.v_j[xlim1:xlim2],norm_flag)
    diff[45-1]=norm(v1.I_j[xlim1:xlim2]-v2.I_j[xlim1:xlim2],norm_flag)/norm(v1.I_j[xlim1:xlim2],norm_flag)
    diff[46-1]=norm(v1.eNOS[xlim1:xlim2]-v2.eNOS[xlim1:xlim2],norm_flag)/norm(v1.eNOS[xlim1:xlim2],norm_flag)
    diff[47-1]=norm(v1.NO_j[xlim1:xlim2]-v2.NO_j[xlim1:xlim2],norm_flag)/norm(v1.NO_j[xlim1:xlim2],norm_flag)
    
    # Wall Mechanics
    diff[48-1]=norm(v1.Mp[xlim1:xlim2]-v2.Mp[xlim1:xlim2],norm_flag)/norm(v1.Mp[xlim1:xlim2],norm_flag)
    diff[49-1]=norm(v1.AMp[xlim1:xlim2]-v2.AMp[xlim1:xlim2],norm_flag)/norm(v1.AMp[xlim1:xlim2],norm_flag)
    diff[50-1]=norm(v1.AM[xlim1:xlim2]-v2.AM[xlim1:xlim2],norm_flag)/norm(v1.AM[xlim1:xlim2],norm_flag)
    diff[51-1]=norm(v1.R[xlim1:xlim2]-v2.R[xlim1:xlim2],norm_flag)/norm(v1.R[xlim1:xlim2],norm_flag)
    
    return diff

def compare_results_a(a1,a2,norm_flag):
    from numpy.linalg import norm
    import import_mat_files as im
    import math
    
    dataset = 'tots_NODRUG_pre'
    fig_hemo_whiskerOptoComparison, LNAME_HET_Data = im.import_Berwick_HET_LNAME_Data(dataset, area='Whisker')

    
    diff_HBO=norm(a1.HBO_N-a2.HBO_N,norm_flag)/norm(a1.HBO_N,norm_flag)
    diff_HBR=norm(a1.HBR_N-a2.HBR_N,norm_flag)/norm(a1.HBR_N,norm_flag)
    
    
    return [diff_HBO, diff_HBR]