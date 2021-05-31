%% Generalize Gear’s method for computing the flow of a viscoelastic fluid %%
%% Author: P. Donald Ariel %%

ninf=10;                                            %% n_infinity
h=0.02;                                             %% step size
n=(ninf)/h+3;                                       %% number of iterations
eta=-2*h:h:ninf;                                    %% Mesh Points
k=[0 0.01 0.05 0.1 0.2 0.3];                %% viscoelastic fluid parameter
s=[1.232587 1.2069 1.1067 0.9939 0.8171 0.6977];    %% f''(0)=Q(0)=s
c(1)=0;                                             %% Constant c0
c(2)=0;                                             %% Constant c1
epsilon=10^(-3);                                    %% error tolerance
err=1;                                              %% current error
Z=zeros(4,n);
for y=1:5
    %while err>epsilon
        f=zeros(1,n);             %% function f
        fp=zeros(1,n);            %% function fp (predicted value of f)
        P=zeros(1,n);             %% function f'=P
        Pp=zeros(1,n);            %% function fp (predicted value of P)
        Q=zeros(1,n);             %% function f''=Q
        Qp=zeros(1,n);            %% function fp (predicted value of Q)
        c(3)=s(y)/2;                         %% constant c2
        c(4)=-1*(1+4*k(y)*(c(3)^2))/6;       %% constant c3
        c(5)=0;                              %% constant c4
        c(6)=(c(3)^2)/30;                    %% constant c5
    
        %%%%----------Function values at first 5 mesh points------------%%%

        for i=1:6
            for j=1:6
                    f(i)=f(i)+c(j)*(((i-3)*h)^(j-1));
            end
        end
        for i=1:5
            for j=1:5
                P(i)=P(i)+j*c(j+1)*(((i-3)*h)^(j-1));
            end
        end
        for i=1:5
            for j=1:4
                Q(i)=Q(i)+j*(j+1)*c(j+2)*(((i-3)*h)^(j-1));
            end
        end

        %%%%------------Fourth-order Predictor-Corrector----------------%%%

        for i=5:n-1
            fp(i+1)=(-10*f(i)+18*f(i-1)-6*f(i-2)+f(i-3)+12*h*P(i))/3;
            Pp(i+1)=(-10*P(i)+18*P(i-1)-6*P(i-2)+P(i-3)+12*h*Q(i))/3;
            Qp(i+1)=(k(y)*f(i)*(15*Q(i)+4*Q(i-1)-14*Q(i-2)+6*Q(i-3)-Q(i-4))+h*(1-2*k(y)*P(i))*(-10*Q(i)+18*Q(i-1)-6*Q(i-2)+Q(i-3))-12*(h^2)*(f(i)*Q(i)+1-(P(i)^2)+k(y)*(Q(i)^2)))/(10*k(y)*f(i)+3*h*(1-2*k(y)*P(i)));
            f(i+1)=(48*f(i)-36*f(i-1)+16*f(i-2)-3*f(i-3)+12*h*Pp(i+1))/25;
            P(i+1)=(48*P(i)-36*P(i-1)+16*P(i-2)-3*P(i-3)+12*h*Qp(i+1))/25;
            Q(i+1)=(k(y)*fp(i+1)*(154*Q(i)-214*Q(i-1)+156*Q(i-2)-61*Q(i-3)+10*Q(i-4))+h*(1-2*k(y)*Pp(i-1))*(48*Q(i)-36*Q(i-1)+16*Q(i-2)-3*Q(i-3))-12*(h^2)*(fp(i+1)*Qp(i+1)+1-(Pp(i+1)^2)+k(y)*(Qp(i+1)^2)))/(45*k(y)*fp(i+1)+25*h*(1-2*k(y)*Pp(i+1)));
        end

        %%%%------Orthogonal collocation using Legendre Polynomial------%%%

        D=optimvar('D',4);  %% Lengendre Constants C
        xi=exp(-1*eta);  %% From serth's paper the transform from eta to xi
        %syms f(xi);
        fun=xi-1-log(xi)+power(xi-1,2).*(D(1)+D(2)*(2.*xi-1)+D(3)*(6*power(xi,2)-6.*xi+1)+D(4)*(20*power(xi,3)-30*power(xi,2)+12.*xi-1));
        %fun=-xi.*(1-(1./xi)+2*(xi-1).*(D(1)+D(2)*(2.*xi-1)+D(3)*(6*power(xi,2)
        %-6.*xi+1)+D(4)*(20*power(xi,3)-30*power(xi,2)+12.*xi-1))+
        %power(xi-1,2).*(2*D(2)+D(3)*(12.*xi-6)+D(4)*(60*power(xi,2)-60.*xi+12)))-1;
        obj=sum((fun-f).^2);    %% Objective function
        lsqproblem=optimproblem("Objective",obj); 
        x0.D=[3/2 1/2 3/2 1/2];
        [sol,fval] = solve(lsqproblem,x0);
        responsedata = evaluate(fun,sol);
        show(lsqproblem)
        fns = fieldnames(sol);
        %A.(fns{3})
        g=sol.(fns{1});
        s(y)=(1+2*(g(1)+g(2)+g(3)+g(4))); %% Updated value of the constant
        %g=diff(diff(responsedata));
        %plot(g);
        %s=g;
        err=evaluate(obj,sol);
        x=evaluate(fun,sol);                %% Updated error
        %figure(1)
        %plot(x)
        %hold on 
        %plot(f)
        %plot(P,f,'r*',f,responsedata,'b-')
        %legend('Original Data','Fitted Curve')
        %xlabel 't'
        %ylabel 'Response'
        %title("Fitted Response")%%
    %end
    Z(y,:)=P;
end
figure(1)
title("Variation of f' with eta for various values of k. Flow near a stagnation point") ;
xlabel("eta") ;
ylabel("f'(eta)") ;
for y=1:5
    plot(eta,Z(y,:));
    hold on;
end
%figure(2)
%plot(eta,f);
%figure(3)
%plot(eta,P);
%figure(4)
%plot(eta,Q);