%% Evaluating the unsteady MHD micropolar fluid flow past 
%% stretching/shirking sheet with heat source and thermal radiation: 
%% Implementing fourth order predictor–corrector FDM
%% Author: M.M. Khader and R.P. Sharma %%

%%%---------------------------Basic Constants---------------------------%%%

Sc=0.2;                   %% Schmidt number
tau=1;                     %% unsteadiness parameter
delta=3;                   %% Thermal Buoyancy parameter
M=1;                       %% Concentration buoyancy parameter
R=0;                       %% Thermal radiation parameter
Ec=0.75;                      %% Eckart Number
K=1;                       %% Micropolar parameter
A1=0.1;                      %% A* 
B1=0.1;                      %% B*
delta1=0;                  %% delta*
Pr=0.71;                    %% Prandtl Number
xinf=12;                    %% n_infinity
hstep=0.02;                %% step size
N=(xinf)/hstep+3;          %% number of iterations
eps=-2*hstep:hstep:xinf;   %% Mesh Points
S=0.5;                      %% Suction/injection paramter  
lambda=1;                %% stretching/shrinking parameter
l1=zeros(1,N);             %% f''(0)=l1
l2=zeros(1,N);             %% h'(0)=l2
l3=zeros(1,N);             %% theta'(0)=l3
l4=zeros(1,N);             %% phi'(0)=l4

%%%---------------------------Initial Constants-------------------------%%%

l1(1)=-0.482185;
l2(1)=-1;
l3(1)=-0.410466;
l4(1)=-0.533872;
%l1(2)=-1.1;
epsilon=10^(-5);
err=1;               %error for the values of l1(k), l2(k), l3(k) and l4(k)
k=1;
%%%-----------------Equations for constant determining------------------%%%
while err>epsilon && k<2
    a(1)=S;                        %% a0
    a(2)=lambda;                   %% a1
    a(3)=l1(k)/2;                  %% a2
    b(1)=-1*l1(k)/2;               %% b0
    b(2)=l2(k);                    %% b1  
    c(1)=1;                        %% c0
    c(2)=l3(k);                    %% c1     
    d(1)=1;                        %% d0
    d(2)=l4(k);                    %% d1
    a(4)=(M*lambda-delta1-delta-K*l2(k)+tau*lambda+lambda^2-S*l1(k))/(6*(1+K));      
    b(3)=(2*K*a(3)+2*K*b(1)+3*tau*b(1)/2+b(1)*a(2)-a(1)*b(2))/(2*(1+K/2));           
    a(5)=(2*M*a(3)-delta1*d(2)-delta*c(2)-2*K*b(3)+3*tau*a(3)+2*a(2)*a(3)-6*a(1)*a(4))/(24*(1+K));       
    b(4)=(6*K*a(4)+2*K*b(2)+2*tau*b(2)+2*a(3)*b(1)-2*a(1)*b(3))/(6*(1+K/2));                              
    c(3)=(2*Pr*tau*c(1)-Pr*(a(1)*c(2)-a(2)*c(1))-4*Pr*Ec*(1+K)*a(3)*a(3)-Pr*M*Ec*a(2)*a(2)-A1*a(2)-B1*c(1))/(2*(1+4*R/3));            
    c(4)=(5*tau*Pr*c(2)/2-24*Pr*Ec*(1+K)*a(3)*a(4)-4*Pr*Ec*M*a(2)*a(3)-2*A1*a(3)-B1*c(2)-2*Pr*(a(1)*c(3)-a(3)*c(1)))/(6*(1+4*R/3)); 
    c(5)=(6*Pr*tau*c(3)-2*Pr*Ec*(1+K)*(48*a(3)*a(5)+36*a(4)*a(4))-Pr*(6*a(1)*c(4)-6*a(4)*c(1))-4*Pr*M*Ec*(2*a(3)*a(3)+3*a(2)*a(3))-6*A1*a(4)-2*B1*c(3))/(24*(1+4*R/3));  
    d(3)=(Sc*tau*d(1)-Sc*(a(1)*d(2)-d(1)*a(2)))/2;                                                         
    a(6)=(6*M*a(4)-2*(delta*c(3)+delta1*d(3))-6*K*b(4)+12*a(4)*tau+4*a(3)*a(3)-24*a(1)*a(5))/(60*(1+K)); 
    c(6)=(21*Pr*tau*c(4)-Pr*(12*a(2)*c(4)-12*a(4)*c(2)+24*c(1)*c(5)-24*a(5)*c(1))-2*Pr*Ec*(1+K)*(120*a(3)*a(6)+432*a(4)*a(5))-2*Pr*M*Ec*(36*a(3)*a(4)+24*a(2)*a(5))-24*A1*a(5)-6*B1*c(4))/(60*(1+4*R/3)); 
    d(4)=(3*Sc*tau*d(2)/2-2*Sc*(a(1)*d(3)-d(1)*a(3)))/6;                                                   
    d(5)=(4*Sc*tau*d(3)-Sc*(2*a(2)*d(3)+6*a(1)*d(4)-6*d(1)*a(4)-2*d(2)*a(3)))/24;                          
    d(6)=(15*Sc*tau*d(4)-Sc*(12*a(2)*d(4)-12*d(2)*a(4)+24*a(1)*d(5)-24*d(1)*a(5)))/60;                    
    b(5)=(K*(4*b(3)+24*a(5))+5*tau*b(3)+6*a(4)*b(1)+2*a(3)*b(2)-6*a(1)*b(4)-2*a(2)*b(3))/(24*(1+K/2));    
    b(6)=(K*(12*b(4)+60*a(6))+18*b(4)*tau+24*a(5)*b(1)-24*a(1)*b(5)+12*a(4)*b(2)-12*a(2)*b(4))/(60*(1+K/2)); 
    
    %%%--------------------------Matrix Creation------------------------%%%
    
    f=zeros(1,N);
    h=zeros(1,N);
    theta=zeros(1,N);
    phi=zeros(1,N);
    F=zeros(1,N);
    H=zeros(1,N);
    Theta=zeros(1,N);
    Phi=zeros(1,N);
    fs=zeros(1,N);
    hs=zeros(1,N);
    thetas=zeros(1,N);
    phis=zeros(1,N);
    Fs=zeros(1,N);
    Hs=zeros(1,N);
    Thetas=zeros(1,N);
    Phis=zeros(1,N);
    
    %%%%-------------Function values at first 5 mesh points-------------%%%
    
    for i=1:5
        for j=1:6
            f(i)=f(i)+a(j)*(eps(i)^(j-1));
            h(i)=h(i)+b(j)*(eps(i)^(j-1));
            theta(i)=theta(i)+c(j)*(eps(i)^(j-1));
            phi(i)=phi(i)+d(j)*(eps(i)^(j-1));
        end
    end
    for i=1:5
        for j=1:5
            F(i)=F(i)+j*a(j+1)*((eps(i))^(j-1));
            H(i)=H(i)+j*b(j+1)*((eps(i))^(j-1));
            Theta(i)=Theta(i)+j*c(j+1)*((eps(i))^(j-1));
            Phi(i)=Phi(i)+j*d(j+1)*((eps(i))^(j-1));
        end
    end  

    %%%%---------------Fourth-order Predictor-Corrector-----------------%%%
    
    for j=5:N-1
        fs(j+1)=(-10*f(j)+18*f(j-1)-6*f(j-2)+f(j-3)+12*hstep*F(j))/3;
        hs(j+1)=(-10*h(j)+18*h(j-1)-6*h(j-2)+h(j-3)+12*hstep*H(j))/3;
        thetas(j+1)=(-10*theta(j)+18*theta(j-1)-6*theta(j-2)+theta(j-3)+12*hstep*Theta(j))/3;
        phis(j+1)=(-10*phi(j)+18*phi(j-1)-6*phi(j-2)+phi(j-3)+12*hstep*Phi(j))/3;
        Fs(j+1)=(-1*(-30*(1+K)+20*f(j)*hstep-10*hstep*tau*eps(j+1)-24*(hstep^2)*(tau+M))*F(j)+(8*(1+K)+36*f(j)*hstep-18*hstep*tau*eps(j+1))*F(j-1)-(28*(1+K)+12*f(j)*hstep-6*hstep*tau*eps(j+1))*F(j-2)+(12*(1+K)+2*f(j)*hstep-hstep*tau*eps(j+1))*F(j-3)+2*(1+K)*F(j-4)+24*(hstep^2)*(F(j)^2-K*H(j)-delta*theta(j)-delta1*phi(j)))/(20*(1+K)+6*hstep*f(j)-3*hstep*tau*eps(j+1));
        F(j+1)=(-1*(30*(1+K)-96*f(j)*hstep-48*hstep*tau*eps(j+1)-24*(hstep^2)*(tau+M))*F(j)+(8*(1+K)-72*f(j)*hstep+36*hstep*tau*eps(j+1))*F(j-1)-(28*(1+K)-32*f(j)*hstep+16*hstep*tau*eps(j+1))*F(j-2)+(12*(1+K)-6*f(j)*hstep+3*hstep+tau*eps(j+1))*F(j-3)+2*(1+K)*F(j-4)+24*(hstep^2)*((F(j)^2)-K*H(j)-delta*theta(j)-delta1*phi(j)))/(20*(1+K)+50*hstep*f(j)-hstep*tau*eps(j+1));
        Hs(j+1)=(-1*(10+12*hstep*(f(j)-0.5*tau*eps(j+1))/(1+0.5*K))*H(j)+18*H(j-1)-6*H(j-2)+H(j-3)+(K*(3*F(j+1)+10*F(j)-18*F(j-1)+6*F(j-2)-F(j-3))+12*hstep*(F(j)+1.5*tau+2*K)*h(j))/(1+0.5*K))/3;
        Thetas(j+1)=(-1*(10*hstep+12*(hstep^2)*Pr*(f(j)-0.5*tau*(eps(j+1)))/(1+4*R/3))*Theta(j)+18*hstep*Theta(j-1)-6*hstep*Theta(j-2)+hstep*Theta(j-3)-12*Pr*Ec*(1+K)*((F(j+1)/4+5*F(j)/6-3*F(j-1)/2+F(j-2)/2-F(j-3)/12)^2)/(1+4*R/3)+12*(hstep^2)*(Pr*theta(j)*(F(j)+tau)-Pr*M*Ec*(F(j)^2)-A1*F(j)-B1*theta(j))/(1+4*R/3))/(3*hstep);
        Phis(j+1)=(-1*(10-6*hstep*tau*Sc*(eps(j+1))+12*hstep*Sc*f(j))*Phi(j)+18*Phi(j-1)-6*Phi(j-2)+Phi(j-3)+12*hstep*phi(j)*Sc*(tau+F(j)))/3;
        f(j+1)=(48*f(j)-36*f(j-1)+16*f(j-2)-3*f(j-3)+12*hstep*Fs(j+1))/25;
        h(j+1)=(48*h(j)-36*h(j-1)+16*h(j-2)-3*h(j-3)+12*hstep*Hs(j+1))/25;
        theta(j+1)=(48*theta(j)-36*theta(j-1)+16*theta(j-2)-3*theta(j-3)+12*hstep*Thetas(j+1))/25;
        phi(j+1)=(48*phi(j)-36*phi(j-1)+16*phi(j-2)-3*phi(j-3)+12*hstep*Phis(j+1))/25;
        H(j+1)=((48-12*hstep*(f(j)-0.5*tau*(eps(j+1))/(1+0.5*K)))*H(j)-36*H(j-1)+16*H(j-2)-3*H(j-3)+(K*(25*F(j+1)-48*F(j)+36*F(j-1)-16*F(j-2)+3*F(j-3))+12*hstep*(F(j)+1.5*tau+2*K)*h(j))/(1+0.5*K))/25;
        Theta(j+1)=((48*hstep-12*(hstep^2)*Pr*(f(j)-0.5*tau*hstep*(j+2))/(1+4*R/3))*Theta(j)-36*hstep*Theta(j-1)+16*hstep*Theta(j-2)-3*hstep*Theta(j-3)-12*Pr*Ec*(1+K)*((25*F(j+1)/12-4*F(j)+3*F(j-1)-4*F(j-2)/3+F(j-3)/4)^2)/(1+4*R/3)+12*(hstep^2)*(Pr*theta(j)*(F(j)+tau)-Pr*M*Ec*(F(j)^2)-A1*F(j)-B1*theta(j))/(1+4*R/3))/(25*hstep);
        Phi(j+1)=((48+12*hstep*tau*Sc*hstep*(j+2)-12*hstep*Sc*f(j))*Phi(j)-36*Phi(j-1)+16*Phi(j-2)-3*Phi(j-3)+12*hstep*phi(j)*Sc*(tau+F(j)))/25;
    end
    %if(k==1)
    %    F1=F(N);
    %end
    %if(k~=1)
    %    l1(k+1)=l1(k)-(l1(k)-l1(k-1))*F(N)/(F(N)-F1);
    %    err=abs(F(N));
    %    F1=F(N);
    %end
    
    %%%-------------------Serth's method----------------------%%%
    
    D=optimvar('D',4);   %% Lengendre Constants C
    xi=exp(-1*eps);      %% From serth's paper the transformation from eta to xi
    %syms f(xi);
    %fun=S+xi-1+log(xi)+power(log(xi),2).*(D(1)+D(2)*(2.*xi-1)+D(3)*(6*power(xi,2)-6.*xi+1)+D(4)*(20*power(xi,3)-30*power(xi,2)+12.*xi-1));
    %fun=xi.*(1-(1./xi)+2*(xi-1).*(D(1)+D(2)*(2.*xi-1)+D(3)*(6*power(xi,2)-6.*xi+1)+D(4)*(20*power(xi,3)-30*power(xi,2)+12.*xi-1))+power(xi-1,2).*(2*D(2)+D(3)*(12.*xi-6)+D(4)*(60*power(xi,2)-60.*xi+12)));
    fun=S.*(D(1)+D(2)*(2.*xi-1)+D(3)*(6*power(xi,2)-6.*xi+1)+D(4)*(20*power(xi,3)-30*power(xi,2)+12.*xi-1));
    obj=sum((fun-f).^2);    %% Objective function
    lsqproblem=optimproblem("Objective",obj); 
    x0.D=[3/2 1/2 3/2 1/2];
    [sol,fval] = solve(lsqproblem,x0);
    %err=evaluate(obj,sol);
    err=abs(f(N));
    x=evaluate(fun,sol);                  %% Updated error
    %figure(1)
    %plot(x);
    %hold on 
    %plot(f);
    %responsedata = evaluate(fun,sol);
    %show(lsqproblem)
    fns = fieldnames(sol);
    %A.(fns{3})
    g=sol.(fns{1});
    syms y;
    %func=S+y-1+log(y)+power(log(y),2).*(g(1)+g(2)*(2.*y-1)+g(3)*(6*power(y,2)-6.*y+1)+g(4)*(20*power(y,3)-30*power(y,2)+12.*y-1));
    func=S.*(g(1)+g(2)*(2.*y-1)+g(3)*(6*power(y,2)-6.*y+1)+g(4)*(20*power(y,3)-30*power(y,2)+12.*y-1));
    l(y)=diff(func);
    m(y)=diff(diff(func));
    %a=double(((l(exp(-xinf)))/m(exp(-xinf))));
    l1(k+1)=l1(k)-double(((l(exp(-xinf)))/m(exp(-xinf))));
    err=max([err,abs(h(N)),abs(theta(N)),abs(phi(N))]);
    
    %%%---------------Newton-Rhapson Method-------------------%%%
    
    l2(k+1)=l2(k)-h(N)/H(N);
    l3(k+1)=l3(k)-theta(N)/Theta(N);
    %l4(k+1)=l4(k)-phi(N)/Phi(N);
    l4(k+1)=l4(k);
    k=k+1;
end
Fi=smooth(F);
figure(2);
plot(eps,F);
%plot(eps,Fi,'b');
figure(3)
plot(eps,h);
figure(4)
plot(eps,theta);
figure(5)
plot(eps,phi);