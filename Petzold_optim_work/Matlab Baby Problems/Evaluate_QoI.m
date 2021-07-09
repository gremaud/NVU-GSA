function [QoI,t,x] = Evaluate_QoI(A,x0,tspan,QoI_flag)
    odefun = @(t,x) A*x;
    [t,x] = ode45(odefun,tspan,x0);
    
    switch QoI_flag
        case 1
            QoI=x(end,end);
        case 2
            QoI=x(end,1);
    end
end

