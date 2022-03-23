from single_evaluation import single_eval
import numpy as np
import csv

x1=-5
x2=5

v,a,t = single_eval('Robin_Pre2','pre',0)

x1_marker=np.where(t>=x1)
x1_marker=x1_marker[0][0]

x2_marker=np.where(t>=x2)
x2_marker=x2_marker[0][0]


iteration='optim_in_progress'
data_choice='pre'
change_index=0
h=0.1

X=np.zeros([len(t[x1_marker:x2_marker]),234])

for i in range(234):
    V_cPlus=np.ones(234)
    V_cMinus=np.ones(234)
    
    V_cPlus[i]=V_cPlus[i]+h
    V_cMinus[i]=V_cMinus[i]-h

    with open('./Samples/Temp_optim.csv','w') as csv_file:
        csv_writer=csv.writer(csv_file,delimiter=',')
        for j in range(234):
            csv_writer.writerow(V_cPlus[[j]])
    v,a,t = single_eval(iteration,data_choice,change_index)
    QoIPlus=a.HBO_N[x1_marker:x2_marker]

    with open('./Samples/Temp_optim.csv','w') as csv_file:
        csv_writer=csv.writer(csv_file,delimiter=',')
        for j in range(234):
            csv_writer.writerow(V_cMinus[[j]])
    v,a,t = single_eval(iteration,data_choice,change_index)
    QoIMinus=a.HBO_N[x1_marker:x2_marker]
    
    QoI_diff=QoIPlus-QoIMinus
    QoI_deriv=QoI_diff/(2*h)
    
    X[:,i]=QoI_deriv
    print(i)
    
XT=np.transpose(X)
F=np.matmul(XT,X)