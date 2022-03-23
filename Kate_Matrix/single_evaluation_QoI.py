
def single_eval(QoI_flag,Param_Index,change_value):
    
    iteration=[Param_Index,change_value]
    data_choice='pre'
    
    import sys
    sys.path.append("./Model_Codes/")
    
    # Import modules
    from scipy.integrate import odeint
    from scipy import io
    from scipy.interpolate import interp1d
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
    import import_mat_files_no_fig as im
    
    from multipliers import V_IC
    
    if data_choice == 'pre':
        p=p_function(iteration,'normal')
    elif data_choice == 'post':
        p=p_function(iteration,'LNAME')
    
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
    u, solver_details = odeint(func, u0, t, args=(p,idx,input_data,total_waveform,start_time), hmax=1e2, full_output=1)
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
    a = set_algebraic_variables(p,u, t, idx, input_data, total_waveform, 'out function')
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
    
    if clean_flag != 1:
        QoI=100
    else:
        

        # Set custom xlims if wanted
        if data_choice == 'pre': 
            dataset = 'tots_LNAME_pre'
        elif data_choice == 'post':
            dataset = 'tots_LNAME_post'
            

        #sys.path.remove("./Model_Codes/")

        
    

        Data = im.import_Berwick_HET_LNAME_Data(dataset, area='Whisker')
        
        if QoI_flag==0:
            interpolator=interp1d(time, a.HBO_N)
            HBO_interp=interpolator(Data.time)
            Error_HBO=HBO_interp-Data.HbOwhisk_mean
            Error_HBO=list(map(lambda x:x*x,Error_HBO))
            Error=Error_HBO
        elif QoI_flag==1:
            interpolator=interp1d(time, a.HBR_N)
            HBR_interp=interpolator(Data.time)
            Error_HBR=HBR_interp-Data.HbRwhisk_mean
            Error_HBR=list(map(lambda x:x*x,Error_HBR))
            Error=Error_HBR
    
    
        QoI=Error

    
    
    return QoI