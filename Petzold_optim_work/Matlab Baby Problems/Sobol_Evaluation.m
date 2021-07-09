function QoI = Sobol_Evaluation(sobol_vector,mirror_flag,inside_parameters,outside_parameters,A,x0,tspan,QoI_flag)
    n_inside=length(inside_parameters);
    n_outside=length(outside_parameters);
    A_temp=A;
    for i=1:n_inside
        current_param=inside_parameters{i};
        A_temp(current_param(1),current_param(2))=A_temp(current_param(1),current_param(2))*sobol_vector(1,i);
        if mirror_flag==true
            A_temp(current_param(2),current_param(1))=A_temp(current_param(2),current_param(1))*sobol_vector(1,i);
        else
            A_temp(current_param(2),current_param(1))=A_temp(current_param(2),current_param(1))*sobol_vector(1,n_inside+n_outside+i);
        end
    end
    for i=1:n_outside
        current_param=outside_parameters{i};
        A_temp(current_param(1),current_param(2))=A_temp(current_param(1),current_param(2))*sobol_vector(1,n_inside+i);
        if mirror_flag==true
            A_temp(current_param(2),current_param(1))=A_temp(current_param(2),current_param(1))*sobol_vector(1,n_inside+i);
        else
            A_temp(current_param(2),current_param(1))=A_temp(current_param(2),current_param(1))*sobol_vector(1,2*n_inside+n_outside+i);
        end
    end
    
    QoI = Evaluate_QoI(A_temp,x0,tspan,QoI_flag);
end

