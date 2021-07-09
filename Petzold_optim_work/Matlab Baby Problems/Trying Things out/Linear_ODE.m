clear
close all

N=25; %number of ODEs/State Variables
Int_max=5; %largest value to appear in the ODE equations
Int_min=-5; %smallest value to appear in the ODE equations
Max_population=10; %Sum of initial conditions
tf=10; %final time

A=randi([Int_min,Int_max],N,N);
for i=1:N
    A(i,i)=0;
    for j=i+1:N
        A(j,i)=-A(i,j);
    end
end

x0=rand(1,N);
x0=x0*Max_population/norm(x0,1);

odefun = @(t,x) A*x;
[t,x] = ode45(odefun,[0,tf],x0);

figure(1)
plot(t,x)

figure(2)
plot(t,x)
xlim([tf-1,tf])

figure(3)
spy(A)