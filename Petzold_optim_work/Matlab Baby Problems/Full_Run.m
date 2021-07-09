% Start with a clean environment
close all
clear

%Initialize varriables
N_list=[4,4,4];
Int_min=-5;
Int_max=5;
Num_connections=1;
Max_population=1;
tspan=0:.01:2;
QoI_flag = 1;   %1-> last value of last x
                %2-> last value of first x
Num_parameters_inside=0;
Num_parameters_outside=Num_connections;

%% Create Original System
tic
[A,x0,list_of_inside,list_of_outside] = Generate_system(N_list,Int_min,Int_max,Num_connections,Max_population);
figure(1)
spy(A)
[nominal_QoI,nominal_t,nominal_x] = Evaluate_QoI(A,x0,tspan,QoI_flag);

inside_parameters=cell(Num_parameters_inside,1);
outside_parameters=cell(Num_parameters_outside,1);

for i=1:Num_parameters_inside
   index=randi([1,length(list_of_inside)]);
   temp=list_of_inside(index);
   list_of_inside(index)=[];
   temp=temp{1};
   new_inside=[temp(2)+sum(N_list(1:temp(1)-1)), temp(3)+sum(N_list(1:temp(1)-1))];
   inside_parameters{i}=new_inside;
end
for i=1:Num_parameters_outside
   index=randi([1,length(list_of_outside)]);
   temp=list_of_outside(index);
   list_of_outside(index)=[];
   temp=temp{1};
   new_outside=temp;
   outside_parameters{i}=new_outside;
end
x=toc;

clear i index Max_population new_inside new_outside temp
disp(['Original system created, time taken: ', num2str(x)])
%% Component Based Analysis
tic

mirror_flag=true; %should (y,x) be paired with (x,y)
component_analysis_flag='QoI';  %'QoI'=Analyze based on just QoI value
                                %'States'=Analyze base on all state values

                                
Component_analysis_values=zeros(Num_parameters_inside+Num_parameters_outside,1);
for i=1:Num_parameters_inside
    current_component=inside_parameters{i};
    A_temp=A;
    A_temp(current_component(1),current_component(2))=0;
    if mirror_flag==true
        A_temp(current_component(2),current_component(1))=0;
    end
    [reaction_QoI,reaction_t,reaction_x] = Evaluate_QoI(A_temp,x0,tspan,QoI_flag);
   
    switch component_analysis_flag
        case 'QoI'
            Component_analysis_values(i)=abs((reaction_QoI-nominal_QoI)/(nominal_QoI));
        case 'States'
            
    end
end
for i=1:Num_parameters_outside
    current_component=outside_parameters{i};
    A_temp=A;
    A_temp(current_component(1),current_component(2))=0;
    if mirror_flag==true
        A_temp(current_component(2),current_component(1))=0;
    end
    [reaction_QoI,reaction_t,reaction_x] = Evaluate_QoI(A_temp,x0,tspan,QoI_flag);
   
    switch component_analysis_flag
        case 'QoI'
            Component_analysis_values(Num_parameters_inside+i)=abs((reaction_QoI-nominal_QoI)/(nominal_QoI));
        case 'States'
            
    end
end

[temp,I]=sort(Component_analysis_values);
parameters_sorted_by_componet=cell(Num_parameters_inside+Num_parameters_outside,2);
for i=1:Num_parameters_inside+Num_parameters_outside
    if I(i)>Num_parameters_inside
        parameters_sorted_by_componet{i,1}=outside_parameters{I(i)-Num_parameters_inside};
    else
        parameters_sorted_by_componet{i,1}=inside_parameters{I(i)};
    end
    parameters_sorted_by_componet{i,2}=temp(i);
end
x=toc;

clear A_temp mirror_flag component_analysis_flag current_component i I ...
    reaction_QoI reaction_t reaction_x temp
disp(['Component Based Analysis Complete, time taken: ', num2str(x)])
%% Sobol Analysis
tic
mirror_flag=true; %should (y,x) be paired with (x,y)
dist_width=0.1;

N_Sobol=1e1;


num_params=Num_parameters_inside+Num_parameters_outside;
if mirror_flag==false
    num_params=num_params*2;
end

sobol_A=zeros(N_Sobol,num_params);
sobol_B=zeros(N_Sobol,num_params);
for i=1:N_Sobol
    sobol_A(i,:)=1-dist_width+2*dist_width*rand(1,num_params);
    sobol_B(i,:)=1-dist_width+2*dist_width*rand(1,num_params);
end

sobol_C_cell=cell(1,num_params);
for i=1:num_params
    sobol_C_temp=sobol_B;
    sobol_C_temp(:,i)=sobol_A(:,i);
    sobol_C_cell{i}=sobol_C_temp;
end

f= @(x) Sobol_Evaluation(x,mirror_flag,inside_parameters,outside_parameters,A,x0,tspan,QoI_flag);

yC=cell(1,num_params);
for i=1:N_Sobol
    yA(i)=f(sobol_A(i,:));
    yB(i)=f(sobol_B(i,:));
    for j=1:num_params
        sobol_C_temp=sobol_C_cell{j};
        yC_temp=yC{j};
        yC_temp(i)=f(sobol_C_temp(i,:));
        yC{j}=yC_temp;
    end
end

f02=N_Sobol^(-2)*sum(yA)*sum(yB);

parfor i=1:num_params
    numS=1/N_Sobol*yA*yC{i}'-f02;
    numST=1/N_Sobol*yB*yC{i}'-f02;
    denom=1/N_Sobol*(yA*yA')-f02;
    S(i)=numS/denom;
    ST(i)=1-(numST/denom);
end

[temp,I]=sort(ST);
parameters_sorted_by_total_sobol=cell(Num_parameters_inside+Num_parameters_outside,2);
for i=1:Num_parameters_inside+Num_parameters_outside
    if I(i)>Num_parameters_inside
        parameters_sorted_by_total_sobol{i,1}=outside_parameters{I(i)-Num_parameters_inside};
    else
        parameters_sorted_by_total_sobol{i,1}=inside_parameters{I(i)};
    end
    parameters_sorted_by_total_sobol{i,2}=temp(i);
end
[temp,I]=sort(S);
parameters_sorted_by_sobol=cell(Num_parameters_inside+Num_parameters_outside,2);
for i=1:Num_parameters_inside+Num_parameters_outside
    if I(i)>Num_parameters_inside
        parameters_sorted_by_sobol{i,1}=outside_parameters{I(i)-Num_parameters_inside};
    else
        parameters_sorted_by_sobol{i,1}=inside_parameters{I(i)};
    end
    parameters_sorted_by_sobol{i,2}=temp(i);
end

x=toc;

clear dist_width f f02 i I j mirror_flag N_Sobol num_params sobol_A ...
    sobol_B sobol_C_cell sobol_C_temp temp yA yB yC yC_temp

disp(['Sobol Based Analysis Complete, time taken: ', num2str(x)])
disp(['Sum of Sobol: ', num2str(sum(S))])

%% Compare the two
parameters_comparison=cell(Num_parameters_inside+Num_parameters_outside,4);
for i=1:Num_parameters_inside+Num_parameters_outside
    if i>Num_parameters_inside
        parameters_comparison{i,1}=outside_parameters{i-Num_parameters_inside};
    else
        parameters_comparison{i,1}=inside_parameters{i};
    end
    for j=1:Num_parameters_inside+Num_parameters_outside
        if isequal(parameters_comparison{i,1},parameters_sorted_by_componet{j,1})
            parameters_comparison{i,2}=parameters_sorted_by_componet{j,2};
        end
        if isequal(parameters_comparison{i,1},parameters_sorted_by_sobol{j,1})
            parameters_comparison{i,3}=parameters_sorted_by_sobol{j,2};
        end
        if isequal(parameters_comparison{i,1},parameters_sorted_by_total_sobol{j,1})
            parameters_comparison{i,4}=parameters_sorted_by_total_sobol{j,2};
        end
    end        
end

