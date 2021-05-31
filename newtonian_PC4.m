
%% Fourth-order predictor-corrector FDM for the effect of viscous dissipation and Joule heating on the Newtonian fluid flow %%
%% Author: M.M. Khader 

Ec=0.1;              %% Local Eckert Number
a1=0.3;              %% a*=temperature-dependent heat source/sink parameter
b1=0.3;              %% b*=space-dependent heat source/sink parameter
R=0.5;               %% R=Radiation paprameter
Pr=1.0;              %% Prandtl number
M=0;                 %% M=local Magnetic parameter
N=50;                %% Number of iterations
h=5/N;               %% Step size
fw=1.5;              %% f=fw at n=0
l=1.621397;          %% theta(0)=l (initial condition)
S=-2.20384;          %% f''(0)=s (initial condition)
x0=4;                %% secant method starting value for root range
x1=5;                %% secant method ending value for root range
epsilon=10^(-4);     %% error tolerance
err=abs(x1-x0);      %% secant method error
while err > epsilon  
    min1=10^(8);
    %for S=-1:0.001:1
    f=zeros(1,N+1);     
    theta=zeros(1,N+1);
    F=zeros(1,N+1);      %% f'= F (velocity profile)
    Theta=zeros(1,N+1);  %% theta'= Theta
    fs=zeros(1,N+1);     %% fs=fp (Predicted values)
    thetas=zeros(1,N+1); %% thetas=thetap (Predicted values)
    Fs=zeros(1,N+1);     %% Fs=Fp (Predicted values)
    Thetas=zeros(1,N+1); %% Thetas=Thetap (Predicted values)
    c(1)=fw;             %% Constant c0
    c(2)=1;              %% Constant c1
    d(1)=l;              %% Constant d0
    d(2)=-1;             %% Constant d1
    c(3)=S/2;            %% Constant c2
    c(4)=(M*c(2)+2*(c(2)^2)-2*c(1)*c(3))/6;                                                
    c(5)=(2*M*c(3)-6*c(1)*c(4)+6*c(2)*c(3))/24;     
    c(6)=(6*M*c(4)-24*c(1)*c(5)+12*c(2)*c(4)+12*(c(3)^2))/120;  
    d(3)=-1*(c(1)*d(2)+c(2)*d(1)+(a1+b1*d(1))/Pr+Ec*(4*(c(3)^2)+
    M*(c(2)^2)))*Pr/(2*(1+R));                         
    d(4)=-1*(2*(c(1)*d(3)+c(3)*d(1)+c(2)*d(2))+(-a1+b1*d(2))/Pr+
    2*Ec*(12*c(3)*c(4)+2*M*c(2)*c(3)))*Pr/(6*(1+R));  
    d(5)=-1*(6*(c(1)*d(4)+c(4)*d(1)+c(3)*d(2)+d(3)*c(2))+(a1+2*b1*d(3))+
    2*Ec*(36*(c(4)^2)+12*c(3)*c(4)+M*(4*(c(3)^2)+6*c(2)*c(4))))*Pr/(24*(1+R));                     
    d(6)=-1*(24*(c(2)*d(4)+c(4)*d(2)+c(3)*d(3)+c(1)*d(5)+d(1)*c(5))+
    (-a1+6*b1*d(4))/Pr+2*Ec*(288*c(4)*c(5)+36*(c(4)^2)+240*c(3)*c(6)+
    M*(36*c(3)*c(4)+24*c(2)*c(5))))*Pr/(120*(1+R));           
    %% Function values at first 5 mesh points %%
    for i=1:5
        for j=1:6
            f(i)=f(i)+c(j)*(((i-3)*h)^(j-1));
            theta(i)=theta(i)+d(j)*(((i-3)*h)^(j-1));
        end
    end
    for i=1:5
        for j=1:5
            F(i)=F(i)+j*c(j+1)*(((i-3)*h)^(j-1));
            Theta(i)=Theta(i)+j*d(j+1)*(((i-3)*h)^(j-1));
        end
    end
    
    %%---------------Fourth-order Predictor-Corrector--------------%%
    for i=5:N
        fs(i+1)=(-10*f(i)+18*f(i-1)-6*f(i-2)+f(i-3)+12*h*F(i))/3;
        thetas(i+1)=(-10*theta(i)+18*theta(i-1)-6*theta(i-2)+theta(i-3)+
        12*h*Theta(i))/3;
        Fs(i+1)=((15-10*h*f(i))*F(i)+(4+18*h*f(i))*F(i-1)-1*(14+6*h*f(i))*F(i-2)+
        (6+h*f(i))*F(i-2)-F(i-4)+24*(h^2)*(F(i)^2)+12*(h^2)*M*F(i))/(10+3*h*f(i));
        Thetas(i+1)=(-1*(10+(12*h*Pr*f(i))/(1+R))*Theta(i)-18*Theta(i-1)+
        6*Theta(i-2)-Theta(i-3)-(12*h/(1+R))*(a1*exp(-h*i)+b1*theta(i))-
        12*h*theta(i)*Theta(i))/3-((4*h*Pr*Ec)/(1+R))*(((F(i+1)/4+5*F(i)/6-
        3*F(i-1)/2+F(i-2)/2-F(i-3)/12)^2)/(h^2)+M*(F(i)^2));
        f(i+1)=(48*f(i)-36*f(i-1)+16*f(i-2)-3*f(i-3)+12*h*Fs(i+1))/25;
        theta(i+1)=(48*theta(i)-36*theta(i-1)+16*theta(i-2)-3*theta(i-3)+
        12*h*Thetas(i+1))/25;
        F(i+1)=(-1*(15+10*h*f(i))*F(i)+(4+18*h*f(i))*F(i-1)-1*(14+6*h*f(i))*F(i-2)+
        (6+h*f(i))*F(i-2)-F(i-4)+24*(h^2)*(F(i)^2)+12*(h^2)*M*F(i))/(10+3*h*f(i));
        Theta(i+1)=(-1*(10+(12*h*Pr*f(i))/(1+R))*Theta(i)-18*Theta(i-1)+6*Theta(i-2)-
        Theta(i-3)-(12*h/(1+R))*(a1*exp(-h*i)+b1*theta(i))-12*h*theta(i)*Theta(i))/3-
        ((4*h*Pr*Ec)/(1+R))*(((F(i+1)/4+5*F(i)/6-3*F(i-1)/2+F(i-2)/2-
        F(i-3)/12)^2)/(h^2)+M*(F(i)^2));
    end
    min1=min(min1,F(N+1));
    if(F(N+1)==min1)
        G=F;
    end
    %%------------------------Secant Method------------------------%%
    if(F(N)<min1)
        min1=F(N);
        ans=F(:,1:N);
    end
    g=@(x)c(2)+2*c(3)*x+3*c(4)*(x^2)+4*c(5)*(x^3)+5*c(6)*(x^4);
    S=(x1*F(x0)-x0*g(x1)/(g(x0)-g(x1));
    x2=(x0*F(round(x1/h))-x1*F(round(x0/h)))/(F(round(x1/h))-F(round(x0/h)));
    err=abs(x1-x0);
    S=(F(round(x2/h))-F(round(x1/h)))/(x2-x1);
end
figure(1)
plot(F);
figure(2)
plot(theta);