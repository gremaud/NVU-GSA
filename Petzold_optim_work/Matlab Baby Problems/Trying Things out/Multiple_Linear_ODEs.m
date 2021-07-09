clear
close all

N_list=[5,5,5]; %number of ODEs/State Variables in each subsystem
N_total=sum(N_list);
Int_max=5; %largest value to appear in the ODE equations
Int_min=-5; %smallest value to appear in the ODE equations
Max_population=10; %Sum of initial conditions
tf=10; %final time
Num_connections=5; %Maximum number of connections between subsystems

A=zeros(N_total,N_total);
for N_iter=1:length(N_list)
    N=N_list(N_iter);
    sub_A=randi([Int_min,Int_max],N,N);
    for i=1:N
        sub_A(i,i)=0;
        for j=i+1:N
            sub_A(j,i)=-sub_A(i,j);
        end
    end
    if N_iter==1
        A(1:N,1:N)=sub_A;
    else
        N_filled=sum(N_list(1:N_iter-1));
        A(N_filled+1:N_filled+N,N_filled+1:N_filled+N)=sub_A;
    end
end

for i=1:Num_connections
    sub_1=randi([1,length(N_list)]);
    if sub_1==length(N_list)
        sub_2=length(N_list);
        sub_1=length(N_list)-1;
    else
        sub_2=sub_1+1;
    end
    
    N1=randi([1,N_list(sub_1)]);
    N2=randi([1,N_list(sub_2)])+sum(N_list(1:sub_2-1));
    value=randi([Int_min,Int_max]);
    while value==0
        value=randi([Int_min,Int_max]);
    end
    A(N1,N2)=value;
    A(N2,N1)=-value;
    
end

% for N_iter=2:length(N_list)
%     N1=sum(N_list(1:N_iter-1));
%     N2=sum(N_list(1:N_iter));
%     value=randi([Int_min,Int_max]);
%     while value==0
%         value=randi([Int_min,Int_max]);
%     end
%     A(N1,N2)=value;
%     A(N2,N1)=-value;
% end

x0=rand(1,N_total);
x0=x0*Max_population/norm(x0,1);


%odefun = @(t,x) A*exp(-t/2)*x;
odefun = @(t,x) A*x;
[t,x] = ode45(odefun,[0,tf],x0);

figure(1)
plot(t,x)

figure(2)
plot(t,x)
xlim([tf-1,tf])

figure(3)
spy(A)