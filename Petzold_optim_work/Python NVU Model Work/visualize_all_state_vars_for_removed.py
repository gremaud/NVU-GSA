import sys
sys.path.append("./Model_Codes/")
from single_evaluation import single_eval
from scipy.interpolate import interp1d
import import_mat_files_no_fig as im
import matplotlib.pyplot as plt
import parameters as p
from matplotlib.patches import Rectangle

from test_switch_functions import compare_results_v
from test_switch_functions import compare_results_a
    
from plotting_functions import get_subplot_dims, get_figure_size
    
import math
import numpy as np
from numpy.linalg import norm


## Tim these three lines are the ones to edit

reactions_to_remove=[4,11,8,52,49,50] # reactions to be removed as a list
xlim1 = -1  #limits /for x values in graph
xlim2 = 20

values='Robin_Pre2'
model='pre'

norm_flag=2
v_a_flag='v'

v_nominal, a_nominal, t =single_eval(values,model,[])
v=[]
a=[]

v_temp, a_temp, t =single_eval(values,model,reactions_to_remove)


if model == 'pre':
    dataset = 'tots_LNAME_pre'
     
    Data = im.import_Berwick_HET_LNAME_Data(dataset, area='Whisker')
    interpolator=interp1d(t, a_temp.HBO_N)
    HBO_interp=interpolator(Data.time)
    Error_HBO=HBO_interp-Data.HbOwhisk_mean
    Error_HBO=sum(map(lambda x:x*x,Error_HBO))

    interpolator=interp1d(t, a_temp.HBR_N)
    HBR_interp=interpolator(Data.time)
    Error_HBR=HBR_interp-Data.HbRwhisk_mean
    Error_HBR=sum(map(lambda x:x*x,Error_HBR))





alphabetical=0
title='All state variables'
figsize=(0,0)
dpi=80
num = 1
show_stimulation=0
plot_list = sorted(vars(v_nominal).keys())



xlim1_idx = np.where(t == xlim1) 
xlim1_idx = int(xlim1_idx[0][0])
 
num_vars = len(plot_list)                               # calculate the number of variables to plot 
row, col = get_subplot_dims(num_vars)                   # get the number of rows and columns for subplot

if figsize == (0,0):
    fig_size = get_figure_size(row,col)                     # get the figure size based on rows and columns
else:
    fig_size = figsize
    
fig = plt.figure(figsize=fig_size, dpi=dpi, edgecolor='k', num=num)   # generate the main figure
plt.tight_layout()
      
i = 1
# loop through and plot the state variables, checking if the variable is in plot_list
loop = sorted(vars(v_nominal).items())  # order alphabetically
loop2= sorted(vars(v_temp).items())
for attr, value in loop:
    if attr in plot_list:
        value2=loop2[i-1][1]
        if num_vars > 1: # if multiple plots then add subplot
            if alphabetical == 1:   # set position according to order of sorted v
                position = i
            else:
                position = plot_list.index(attr) + 1 # find position of variable in plot_list and set the subplot position accordingly (with +1 as indexing begins at zero) 
            ax = fig.add_subplot(row,col,position)
            ax.set_title(attr)
        else:
            ax = fig.gca()
            variable_name = attr
            
        ax.plot(t[xlim1_idx:], value[xlim1_idx:],t[xlim1_idx:], value2[xlim1_idx:])               
        plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)   # stop weird exponential notation
        plt.xlabel('Time [s]')
        plt.xlim(xlim1,xlim2)
        plt.tight_layout()
        if abs((max(value[xlim1_idx:]) - min(value[xlim1_idx:]))) < 1e-6: # check if the array is pretty much constant and change ylim to just below & above the mean if so (otherwise gives a misleading figure)
            plt.ylim( 0.9 * np.mean(value[xlim1_idx:]), 1.1 * np.mean(value[xlim1_idx:]) )
#            if abs(max(value[xlim1_idx:])) < 1e-6 and abs(min(value[xlim1_idx:])) < 1e-6 and abs(np.mean(value[xlim1_idx:])) < 1e-6: # check if the array is zero for all time and if so set ylim to -1,1
#                plt.ylim(-1,1)
        if show_stimulation == 1:
            ymin, ymax = plt.gca().get_ylim()
            plt.gca().add_patch(Rectangle((0, ymin), p.lengthpulse/1e3, ymax-ymin, facecolor='k', alpha=0.1)) # plot an optional grey square where the stimulation is
        i = i + 1
        
# loop through and plot the algebraic variables, checking if the variable is in plot_list  
# loop = sorted(vars(a_nominal).items()) # order alphabetically
# for attr, value in loop:
#     if attr in plot_list:
        
#         if num_vars > 1: # if multiple plots then add subplot
#             if alphabetical == 1:   # set position according to order of sorted a 
#                 position = i
#             else:
#                 position = plot_list.index(attr) + 1 # find position of variable in the plot_list list and set the subplot position accordingly (with +1 as indexing begins at zero)
#             ax = fig.add_subplot(row,col,position)
#             ax.set_title(attr)
#         else:
#             ax = fig.gca()  
#             variable_name = attr
                            
#         if isinstance(value, np.ndarray):   # check if the algebraic variable is an array
#             ax.plot(t[xlim1_idx:], value[xlim1_idx:])
#             if abs((max(value[xlim1_idx:]) - min(value[xlim1_idx:]))) < 1e-6: # check if the array is pretty much constant and change ylim to below & above the mean if so (otherwise gives a misleading figure)
#                 plt.ylim( 0.9 * np.mean(value[xlim1_idx:]), 1.1 * np.mean(value[xlim1_idx:]) )
# #                if abs(max(value[xlim1_idx:])) < 1e-6 and abs(min(value[xlim1_idx:])) < 1e-6 and abs(np.mean(value[xlim1_idx:])) < 1e-6: # check if the array is zero for all time and if so set ylim to -1,1
# #                    plt.ylim(-1,1)
#         else:                               # otherwise it's a constant (not array, must be handled differently)  
#             ax.axhline(y=value)                     
#             if abs(value) < 1e-6:
#                 plt.ylim(-1,1)  # if value is zero then set ylim to -1,1
#             else:
#                 plt.ylim(0.9*value, 1.1*value) # otherwise set ylim to just below & above value
                
#         if show_stimulation == 1:
#             ymin, ymax = plt.gca().get_ylim()
#             plt.gca().add_patch(Rectangle((0, ymin), p.lengthpulse/1e3, ymax-ymin, facecolor='k', alpha=0.1)) # plot an optional grey square where the stimulation is
#         plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
#         plt.xlabel('Time [s]')
#         plt.xlim(xlim1,xlim2)
#         plt.tight_layout()
#         i = i + 1

if not title and num_vars == 1:             # if there is no user specified title and only one variable then make the title the name of the variable
    plt.title(variable_name, y=1)
elif not title and num_vars > 1:            # if there is no user specified title but more than one variable then print with no title
    pass
else:
   fig.suptitle(title, y=1)              # otherwise use the user specified title
   
fig.subplots_adjust(hspace = 1)
plt.show()
temp=fig
name="./State_variables"+''.join(['_'+str(elem) for elem in reactions_to_remove])+'.pdf'
temp.savefig(name, bbox_inches='tight')

            
print('done')


