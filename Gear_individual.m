%% Generalize Gear’s method for computing the flow of a viscoelastic fluid %%
%% Author: P. Donald Ariel %%

ninf=8;          %% n_infinity
h=0.02;            %% step size
n=(ninf)/h+3;      %% number of iterations
eta=-2*h:h:ninf;   %% Mesh Points
k=0;               %% viscoelastic fluid parameter
s=1.2;         %% f''(0)=Q(0)=s
c(1)=0;            %% Constant c0
c(2)=0;            %% Constant c1
epsilon=10^(-2);   %% error tolerance
err=100;             %% current error
while err>epsilon
        f=zeros(1,n);    %% function f
        fp=zeros(1,n);   %% function fp (predicted value of f)
        P=zeros(1,n);    %% function f'=P
        Pp=zeros(1,n);   %% function fp (predicted value of P)
        Q=zeros(1,n);    %% function f''=Q
        Qp=zeros(1,n);   %% function fp (predicted value of Q)
        c(3)=s/2;                            %% constant c2
        c(4)=-1*(1+4*k*(c(3)^2))/6;          %% constant c3
        c(5)=0;                              %% constant c4
        c(6)=(c(3)^2)/30;                    %% constant c5
    
        %%%%-------Function values at first 5 mesh points------------%%%

        for i=1:5
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

        %%%%----------Fourth-order Predictor-Corrector------------%%%

        for i=5:n-1
            fp(i+1)=(-10*f(i)+18*f(i-1)-6*f(i-2)+f(i-3)+12*h*P(i))/3;
            Pp(i+1)=(-10*P(i)+18*P(i-1)-6*P(i-2)+P(i-3)+12*h*Q(i))/3;
            Qp(i+1)=(k*f(i)*(15*Q(i)+4*Q(i-1)-14*Q(i-2)+6*Q(i-3)-Q(i-4))+h*(1-2*k*P(i))*(-10*Q(i)+18*Q(i-1)-6*Q(i-2)+Q(i-3))-12*(h^2)*((f(i)*Q(i))+1-(P(i)^2)+k*(Q(i)^2)))/(10*k*f(i)+3*h*(1-2*k*P(i)));
            f(i+1)=(48*f(i)-36*f(i-1)+16*f(i-2)-3*f(i-3)+12*h*Pp(i+1))/25;
            P(i+1)=(48*P(i)-36*P(i-1)+16*P(i-2)-3*P(i-3)+12*h*Qp(i+1))/25;
            Q(i+1)=(k*fp(i+1)*(154*Q(i)-214*Q(i-1)+156*Q(i-2)-61*Q(i-3)+10*Q(i-4))+h*(1-2*k*Pp(i-1))*(48*Q(i)-36*Q(i-1)+16*Q(i-2)-3*Q(i-3))-12*(h^2)*((fp(i+1)*Qp(i+1))+1-(Pp(i+1)^2)+k*(Qp(i+1)^2)))/(45*k*fp(i+1)+25*h*(1-2*k*Pp(i+1)));
        end
        %figure(2)
        %plot(eta,f);
        %figure(3)
        %plot(eta,P);
        %figure(4)
        %plot(eta,Q);
        %%%%----------Orthogonal collocation using Legendre Polynomial------------%%%

        D=optimvar('D',4);   %% Lengendre Constants C
        xi=exp(-1*eta);      %% From serth's paper the transformation from eta to xi
        %syms f(xi);
        fun=xi-1-log(xi)+power(log(xi),2).*(D(1)+D(2)*(2.*xi-1)+D(3)*(6*power(xi,2)-6.*xi+1)+D(4)*(20*power(xi,3)-30*power(xi,2)+12.*xi-1));
        %fun=xi.*(1-(1./xi)+2*(xi-1).*(D(1)+D(2)*(2.*xi-1)+D(3)*(6*power(xi,2)-6.*xi+1)+D(4)*(20*power(xi,3)-30*power(xi,2)+12.*xi-1))+power(xi-1,2).*(2*D(2)+D(3)*(12.*xi-6)+D(4)*(60*power(xi,2)-60.*xi+12)));
        obj=sum((fun-f).^2);    %% Objective function
        lsqproblem=optimproblem("Objective",obj); 
        x0.D=[3/2 1/2 3/2 1/2];
        [sol,fval] = solve(lsqproblem,x0);
        %err=evaluate(obj,sol);
        err=abs(P(n)-1);
        x=evaluate(fun,sol);                  %% Updated error
        %responsedata = evaluate(fun,sol);
        %show(lsqproblem)
        fns = fieldnames(sol);
        %A.(fns{3})
        g=sol.(fns{1});
        syms y;
        func=y-1-log(y)+power(log(y),2).*(g(1)+g(2)*(2.*y-1)+g(3)*(6*power(y,2)-6.*y+1)+g(4)*(20*power(y,3)-30*power(y,2)+12.*y-1));
        l(y)=diff(func);
        m(y)=diff(diff(func));
        a=double(((l(exp(-ninf))-1)/m(exp(-ninf))));
        s=s-double(((l(exp(-ninf))-1)/m(exp(-ninf))));
        %s=double(power((-m(1)-1)/k,1/2));
        %s=double(l(1));
        %s=(1+2*(g(1)+g(2)+g(3)+g(4)));         %% Updated value of the initial constant
        %g=diff(diff(responsedata));
        %plot(g);
        %s=g;
        %err=max(abs(x-f));
        %figure(1)
        %plot(x)
        %hold on
        %plot(f)
        %figure(5) 
        %fplot(l(y));
        %plot(f)
        %plot(P,f,'r*',f,responsedata,'b-')
        %legend('Original Data','Fitted Curve')
        %xlabel 't'
        %ylabel 'Response'
        %title("Fitted Response")%%
end    
figure(2)
plot(eta,f);
figure(3)
plot(eta,P);
figure(4)
plot(eta,Q);