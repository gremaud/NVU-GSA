import sys
sys.path.append("./Model_Codes/")

#os.chdir('./Model_Codes')

# Import modules
from scipy.integrate import odeint
from scipy import io
import numpy as np
import time
import matplotlib.pyplot as plt
import warnings
from matplotlib.patches import Rectangle

# Local files
from indices import set_indices
from ICs import set_initial_conditions
from ODEsystem import func
from algebraic_variables import set_algebraic_variables
from state_variable_names import initialise_var_names
from normalised_hemo import solve_normalised_hemodynamics
from plotting_functions import plot_variables_singles, plot_variables_samegraph
from model_functions import T_pulses
from parameters import p_function
import import_mat_files as im

from multipliers import V_IC



iteration='Robin_Pre2'
data_choice='pre'
change_index=0



if data_choice == 'pre' :
    p=p_function(iteration,'normal',change_index)
elif data_choice == 'post':
    p=p_function(iteration,'LNAME',change_index)

# Calculate the start time
start_time = time.time()

# Initialise indices, initial conditions and time vector
idx = set_indices()
u0 = set_initial_conditions(idx, V_IC)
t = np.arange(start=0, stop=p.Tend+p.dt, step=p.dt)

# Import the experimental neural input profile from mat file, for InputCase = 'ZhengData', 'ThalamicTrianglesZheng' or 'ZhengFittedParams'
input_data = im.import_Zheng_neural_data(p,t)

# Generate the multiple triangular pulse input profile for the thalamic input T(t), for InputCase = 'ThalamicTriangles' or 'ThalamicTrianglesZheng'
total_waveform = T_pulses(p)

# Solve ODE system
#print("\n Model simulation successfully started\n")
u, solver_details = odeint(func, u0, t, args=(p,idx,input_data,total_waveform,start_time,change_index), hmax=1e2, full_output=1)
if solver_details['message'] == 'Integration successful.':
    del solver_details

# Print total time taken
solve_time = time.time() - start_time
if solve_time < 60.0:
    print("\n\n--- Took {0:.0f} seconds to solve ---".format(solve_time))
else:
    minutes, seconds= divmod(solve_time, 60)
    print("\n\n--- Took {0:.0f} minutes and {1:.0f} seconds to solve ---".format(minutes, seconds))
del start_time, solve_time

# Extract the state variables for easy plotting etc
v = initialise_var_names(u, idx, 'out function')

# Extract the algebraic variables for easy plotting etc
a = set_algebraic_variables(p,u, t, idx, input_data, total_waveform, 'out function',change_index)
del total_waveform, input_data

# Solve for the normalised CBF, CBV, CMRO2, HbO, HbR, HbT and BOLD and put into 'a' (solved separately from odeint as they depend on the steady state conditions)
a = solve_normalised_hemodynamics(p,a,v,t)

# The variables that can be plotted sorted into alphabetical order
state_vars = sorted(vars(v).keys())
#print("\nState variables: \n", state_vars, "\n")
alg_vars = sorted(vars(a).keys())
#print("Algebraic variables: \n", alg_vars,"\n")







if p.startpulse < p.Tend:
    time = t/1e3 - p.startpulse/1e3     # Time in seconds, normalised so t=0 at the start of stimulation
    xlim1 = p.startpulse/1e3 - 10       # Set left xlimit to 10 sec before stimulation begins
else:
    time = t/1e3                        # Time in seconds, not normalised if there is no stimulation
    xlim1 = 0                           # Set left xlimit to 0 if there is no stimulation
    xlim2 = p.Tend/1e3  # Set right xlimit to end of simulation
del t


pulse_marker=np.where(time>=0.0)
pulse_marker=pulse_marker[0][0]

pre_pulse_marker=np.where(time>=-5.0)
pre_pulse_marker=pre_pulse_marker[0][0]

pre_end_marker=np.where(time>=time[-1]-5.0)
pre_end_marker=pre_end_marker[0][0]

radius=v.R
radius = radius/radius[pulse_marker]

clean_flag=0
if not np.all(np.isfinite(radius)):
    clean_flag=-2
elif not np.all(np.isreal(radius)):
    clean_flag=-3
elif np.amax( abs(radius[pre_pulse_marker:pulse_marker+1]-1))>1e-3:
    clean_flag=-4
elif np.amax( abs(radius[pre_end_marker:]-1))>1e-3:
    clean_flag=-5    
elif np.amin(abs(radius[pulse_marker:]))<0.9:
    clean_flag=-6  
else:
    clean_flag=1



# Set custom xlims if wanted
xlim1 = -1
xlim2 = 20
if data_choice == 'pre': 
    dataset = 'tots_NODRUG_pre'
elif data_choice == 'post':
    dataset = 'tots_LNAME_post'
    

fig_hemo_whiskerOptoComparison, LNAME_HET_Data = im.import_Berwick_HET_LNAME_Data(dataset, area='Whisker')

# Plot the hemodynamics for the model, whisker and opto experiments with error bars - 3 subplots for HbO, HbR, HbT
fig_hemo_whiskerOptoComparison = plt.figure(figsize=(14,5), dpi=80, edgecolor='k') 
#
fig_hemo_whiskerOptoComparison.add_subplot(1,3,1)
plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
plt.plot(time, a.HBO_N, label = 'model')
#plt.errorbar(LNAME_HET_Data.time, LNAME_HET_Data.HbOopto_mean, LNAME_HET_Data.HbOopto_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Opto')
plt.errorbar(LNAME_HET_Data.time, LNAME_HET_Data.HbOwhisk_mean, LNAME_HET_Data.HbOwhisk_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Whisk')
plt.title('HbO')
plt.xlabel('Time [s]') 
plt.xlim(-2,20)  
plt.ylim(0.94, 1.1)  
plt.legend(loc='best')
plt.tight_layout()
#
fig_hemo_whiskerOptoComparison.add_subplot(1,3,2)
plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
plt.plot(time, a.HBR_N, label = 'model')
#plt.errorbar(LNAME_HET_Data.time, LNAME_HET_Data.HbRopto_mean, LNAME_HET_Data.HbRopto_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Opto')
plt.errorbar(LNAME_HET_Data.time, LNAME_HET_Data.HbRwhisk_mean, LNAME_HET_Data.HbRwhisk_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Whisk')
plt.title('Reaction '+ str(change_index)+ '; HbR')
plt.xlabel('Time [s]') 
plt.xlim(-2,20)  
plt.ylim(0.94, 1.1)  
plt.legend(loc='best')
plt.tight_layout()
#
fig_hemo_whiskerOptoComparison.add_subplot(1,3,3)
plt.gca().get_yaxis().get_major_formatter().set_useOffset(False) # stop weird exponential notation
plt.plot(time, a.HBT_N, label = 'model')
#plt.errorbar(LNAME_HET_Data.time, LNAME_HET_Data.HbTopto_mean, LNAME_HET_Data.HbTopto_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Opto')
plt.errorbar(LNAME_HET_Data.time, LNAME_HET_Data.HbTwhisk_mean, LNAME_HET_Data.HbTwhisk_std, capsize=3, errorevery=12, elinewidth=1.5, capthick=1, label = 'Whisk')
plt.title('HbT')
plt.xlabel('Time [s]') 
plt.xlim(-2,20)  
plt.ylim(0.94, 1.1)  
plt.legend(loc='best')
plt.savefig('Reaction_figures/Reaction'+str(change_index)+'.pdf')