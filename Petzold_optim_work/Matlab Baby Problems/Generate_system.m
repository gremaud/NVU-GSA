function [A,x0,list_of_inside,list_of_outside] = Generate_system(N_list,Inside_min,Inside_max,Outside_min,Outside_max,Num_connections,Max_population)
    list_of_inside={};
    list_of_outside={};

    N_total=sum(N_list);
    
    A=zeros(N_total,N_total);
    for N_iter=1:length(N_list)
        N=N_list(N_iter);
        sub_A=Inside_min+(Inside_max-Inside_min)*rand(N,N);
        for i=1:N
            sub_A(i,i)=0;
            for j=i+1:N
                sub_A(j,i)=-sub_A(i,j);
                if sub_A(i,j)~=0
                    list_of_inside{end+1}=[N_iter,i,j];
                end
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
            sub_1=1;
        else
            sub_2=sub_1+1;
        end

        N1=randi([1,N_list(sub_1)])+sum(N_list(1:sub_1-1));
        N2=randi([1,N_list(sub_2)])+sum(N_list(1:sub_2-1));
        value=Outside_min+(Outside_max-Outside_min)*rand;
        while value==0
            value=Outside_min+(Outside_max-Outside_min)*rand;
        end
        A(N1,N2)=value;
        A(N2,N1)=-value;
        list_of_outside{end+1}=[N1,N2];

    end
    
    x0=abs(rand(1,N_total));
    x0=x0*Max_population/norm(x0);
    
end

