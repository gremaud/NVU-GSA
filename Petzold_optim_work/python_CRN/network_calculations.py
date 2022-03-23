def QoI(params):
    from scipy.integrate import odeint as solver
    import numpy as np
    
    x0=[500,0,0,50,100,0,0,0,0,0,0]
    
    func=lambda x,t: RHS(x,t,[1,1,1,1,1,1])
    
    t=np.linspace(0, 1, 100)
    sol = solver(func, x0,t)
    
    x0=sol[-1,:]
    func=lambda x,t: RHS(x,t,params)
    t=np.linspace(1, 20, 1000)
    sol = solver(func, x0,t)
    return sol[-1,-1]
    

def RHS(x,t,params):
    M, Mp, Mpp, MAPKK, MKP3, M_MAPKK, Mp_MAPKK, Mpp_MKP3, Mp_MKP3_dep, Mp_MKP3, M_MKP3= x
    k1 =0.02
    k_1= params[0]*1
    k2 =0.01
    k3 =0.032
    k_3= params[1]*1
    k4 =15
    h1 =0.045
    h_1= params[2]*1
    h2 =0.092
    h3 =1
    h_3= params[3]*0.01
    h4 =0.01
    h_4= params[4]*1
    h5 =0.5
    h6 =0.086
    h_6=params[5]*0.001
    
    v1=k1*M*MAPKK-k_1*M_MAPKK
    v2=k2*M_MAPKK
    v3=k3*Mp*MAPKK-k_3*Mp_MAPKK
    v4=k4*Mp_MAPKK
    v5=h1*Mpp*MKP3-h_1*Mpp_MKP3
    v6=h2*Mpp_MKP3
    v7=h3*Mp_MKP3_dep-h_3*Mp*MKP3
    v8=h4*Mp*MKP3-h_4*Mp_MKP3
    v9=h5*Mp_MKP3
    v10=h6*M_MKP3-h_6*M*MKP3
    
    dydt= [v10-v1,v2+v7-v3-v8,v4-v5,v2+v4-v1-v3,v7+v10-v5-v8,v1-v2,v3-v4,v5-v6,v6-v7,v8-v9,v9-v10]
    return dydt


def plot(sol,t):
    import matplotlib.pyplot as plt
    for i in range(10):
        plt.plot(t, sol[:, i])
    plt.show()
    