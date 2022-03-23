import sys
sys.path.append("/Users/rmmorill/Documents/GitHub/NVU-GSA/Petzold_optim_work/Python NVU Model Work/")

from discrete_sobol_with_holes import discrete_sobol_holes_even_probability as sobol
from network_calculations import QoI

import numpy as np
import random
import statistics as stats

list_of_low_reactions=[2, 5, 7, 11, 13, 16]

combos=list()
n = len(list_of_low_reactions)  
for i in range(2**n):
    b = bin(i)[2:]
    l = len(b)
    b = str(0) * (n - l) + b
    combos.append(b)

numpairs=len(combos)
combos_for_ST=combos.copy()

QoIs=np.zeros(numpairs)
list_removed=list()
for i in range(numpairs):
    
    current_removed_index=list(combos[i])
    for j in range(len(current_removed_index)):
        current_removed_index[j]=int(current_removed_index[j])
        
    QoIs[i]=QoI(current_removed_index)


comp_S, comp_ST = sobol(combos,combos,QoIs)


QoI_true=QoI([1,1,1,1,1,1])
QoIs_error=np.zeros(numpairs)

min1_error=QoI_true
min1_index=-1
min2_error=QoI_true
min2_index=-2
for i in range(numpairs):
    
    current_removed_index=list(combos[i])
    for j in range(len(current_removed_index)):
        current_removed_index[j]=int(current_removed_index[j])
    
    QoIs_error[i]= abs(QoI_true-QoIs[i])
    if QoIs_error[i]<min1_error and sum(current_removed_index)<=6-5:
        min1_error=QoIs_error[i]
        min1_index=i
        
    if QoIs_error[i]<min2_error and sum(current_removed_index)<=6-2:
        min2_error=QoIs_error[i]
        min2_index=i



#N=50000
# A=np.zeros([N,6])
# B=np.zeros([N,6])
# for i in range(N):
#     A[i,:]=[random.uniform(0.9,1.1),random.uniform(0.9,1.1),random.uniform(0.9,1.1),random.uniform(0.9,1.1),random.uniform(0.9,1.1),random.uniform(0.9,1.1)]
#     B[i,:]=[random.uniform(0.9,1.1),random.uniform(0.9,1.1),random.uniform(0.9,1.1),random.uniform(0.9,1.1),random.uniform(0.9,1.1),random.uniform(0.9,1.1)]

# C1=np.copy(A)
# C2=np.copy(A)
# C3=np.copy(A)
# C4=np.copy(A)
# C5=np.copy(A)
# C6=np.copy(A)

# C1[:,0]=B[:,0]
# C2[:,1]=B[:,1]
# C3[:,2]=B[:,2]
# C4[:,3]=B[:,3]
# C5[:,4]=B[:,4]
# C6[:,5]=B[:,5]

# fa=np.zeros(N)  
# fb=np.zeros(N)    
# fc1=np.zeros(N)   
# fc2=np.zeros(N)  
# fc3=np.zeros(N)   
# fc4=np.zeros(N)  
# fc5=np.zeros(N)    
# fc6=np.zeros(N)    

# for i in range(N):
#     print(i,' of ', N)
#     fa[i]=QoI(A[i,:])
#     fb[i]=QoI(B[i,:])
#     fc1[i]=QoI(C1[i,:])
#     fc2[i]=QoI(C2[i,:])
#     fc3[i]=QoI(C3[i,:])
#     fc4[i]=QoI(C4[i,:])
#     fc5[i]=QoI(C5[i,:])
#     fc6[i]=QoI(C6[i,:])


# f02=(N**-2)*sum(fa)*sum(fb)
# denom=1/N*np.dot(fa,fa)-f02;

# numST1=1/N*np.dot(fb,fc1)-f02
# numST2=1/N*np.dot(fb,fc2)-f02
# numST3=1/N*np.dot(fb,fc3)-f02
# numST4=1/N*np.dot(fb,fc4)-f02
# numST5=1/N*np.dot(fb,fc5)-f02
# numST6=1/N*np.dot(fb,fc6)-f02

# param_ST=[1-numST1/denom,1-numST2/denom,1-numST3/denom,1-numST4/denom,1-numST5/denom,1-numST6/denom]

print(comp_ST)
print(combos[min1_index])
print(combos[min2_index])
# print(param_ST)