N=25;
f=zeros(1,N+3);
h=zeros(1,N+3);
theta=zeros(1,N+3);
phi=zeros(1,N+3);
F=zeros(1,N+3);
H=zeros(1,N+3);
Theta=zeros(1,N+3);
Phi=zeros(1,N+3);
fs=zeros(1,N);
hs=zeros(1,N);
thetas=zeros(1,N);
phis=zeros(1,N);
Fs=zeros(1,N);
Hs=zeros(1,N);
Thetas=zeros(1,N);
Phis=zeros(1,N);
hstep=6/N;
k=-2:N;
eps=k*hstep;
for i=1:5
    for j=1:6
        f(i)=f(i)+a(j)*(((i-3)*hstep)^(j-1));
        h(i)=h(i)+b(j)*(((i-3)*hstep)^(j-1));
        theta(i)=theta(i)+c(j)*(((i-3)*hstep)^(j-1));
        phi(i)=phi(i)+d(j)*(((i-3)*hstep)^(j-1));
    end
end
for i=1:5
    for j=1:5
        F(i)=F(i)+j*a(j+1)*(((i-3)*hstep)^(j-1));
        H(i)=H(i)+j*b(j+1)*(((i-3)*hstep)^(j-1));
        Theta(i)=Theta(i)+j*c(j+1)*(((i-3)*hstep)^(j-1));
        Phi(i)=Phi(i)+j*d(j+1)*(((i-3)*hstep)^(j-1));
    end
end  

for j=5:N+2
    %F(j)=(f(j+1)/4+(5*f(j))/6-(3*f(j-1))/2+f(j-2)/2-f(j-3)/12)/hstep;
    %H(j)=(h(j+1)/4+(5*h(j))/6-(3*h(j-1))/2+h(j-2)/2-h(j-3)/12)/hstep;
    %Theta(j)=(theta(j+1)/4+(5*theta(j))/6-(3*theta(j-1))/2+theta(j-2)/2-theta(j-3)/12)/hstep;
    %Phi(j)=(phi(j+1)/4+(5*phi(j))/6-(3*phi(j-1))/2+phi(j-2)/2-phi(j-3)/12)/hstep;
    fs(j+1)=(-10*f(j)+18*f(j-1)-6*f(j-2)+f(j-3)+12*hstep*F(j))/3;
    hs(j+1)=(-10*h(j)+18*h(j-1)-6*h(j-2)+h(j-3)+12*hstep*H(j))/3;
    thetas(j+1)=(-10*theta(j)+18*theta(j-1)-6*theta(j-2)+theta(j-3)+12*hstep*Theta(j))/3;
    phis(j+1)=(-10*phi(j)+18*phi(j-1)-6*phi(j-2)+phi(j-3)+12*hstep*Phi(j))/3;
    Fs(j+1)=(-1*(30*(1+K)+20*f(j)*hstep-10*hstep*tau*((j+2)*hstep)-24*(hstep^2)*(tau+M))*F(j)+(8*(1+K)+36*f(j)*hstep-18*hstep*tau*((j+2)*hstep))*F(j-1)-(28*(1+K)+12*f(j)*hstep-6*hstep*tau*((j+2)*hstep))*F(j-2)+(12*(1+K)+2*f(j)*hstep-hstep*tau*((j+2)*hstep))*F(j-3)+2*(1+K)*F(j-4)+24*(hstep^2)*(F(j)^2-K*H(j)-delta*theta(j)-delta1*phi(j)))/(20*(1+K)+6*hstep*f(j)-3*hstep*tau*((j+2)*hstep));
    F(j+1)=(-1*(30*(1+K)-96*f(j)*hstep-48*hstep*tau*(j+2)*hstep-24*(hstep^2)*(tau+M))*F(j)+(8*(1+K)-72*f(j)*hstep+36*hstep*tau*(j+2)*hstep)*F(j-1)-(28*(1+K)-32*f(j)*hstep+16*hstep*tau*(j+2)*hstep)*F(j-2)+(12*(1+K)-6*f(j)*hstep+3*hstep+tau*(j+2)*hstep)*F(j-3)+2*(1+K)*F(j-4)+24*(hstep^2)*((F(j)^2)-K*H(j)-delta*theta(j)-delta1*phi(j)))/(20*(1+K)+50*hstep*f(j)-hstep*tau*(j+2)*hstep);
    Hs(j+1)=(-1*(10+12*hstep*(f(j)-0.5*tau*((j+2)*hstep))/(1+0.5*K))*H(j)+18*H(j-1)-6*H(j-2)+H(j-3)+(K*(3*F(j+1)+10*F(j)-18*F(j-1)+6*F(j-2)-F(j-3))+12*hstep*(F(j)+1.5*tau+2*K)*h(j))/(1+0.5*K))/3;
    Thetas(j+1)=(-1*(10*hstep+12*(hstep^2)*Pr*(f(j)-0.5*tau*((j+2)*hstep))/(1+4*R/3))*Theta(j)+18*hstep*Theta(j-1)-6*hstep*Theta(j-2)+hstep*Theta(j-3)-12*Pr*Ec*(1+K)*((F(j+1)/4+5*F(j)/6-3*F(j-1)/2+F(j-2)/2-F(j-3)/12)^2)/(1+4*R/3)+12*(hstep^2)*(Pr*theta(j)*(F(j)+tau)-Pr*M*Ec*(F(j)^2)-A1*F(j)-B1*theta(j))/(1+4*R/3))/(3*hstep);
    Phis(j+1)=(-1*(10-6*hstep*tau*Sc*((j+2)*hstep)+12*hstep*Sc*f(j))*Phi(j)+18*Phi(j-1)-6*Phi(j-2)+Phi(j-3)+12*hstep*phi(j)*Sc*(tau+F(j)))/3;
    f(j+1)=(48*f(j)-36*f(j-1)+16*f(j-2)-3*f(j-3)+12*hstep*Fs(j+1))/25;
    h(j+1)=(48*h(j)-36*h(j-1)+16*h(j-2)-3*h(j-3)+12*hstep*Hs(j+1))/25;
    theta(j+1)=(48*theta(j)-36*theta(j-1)+16*theta(j-2)-3*theta(j-3)+12*hstep*Thetas(j+1))/25;
    phi(j+1)=(48*phi(j)-36*phi(j-1)+16*phi(j-2)-3*phi(j-3)+12*hstep*Phis(j+1))/25;
    H(j+1)=((48-12*hstep*(f(j)-0.5*tau*((j+2)*hstep)/(1+0.5*K)))*H(j)-36*H(j-1)+16*H(j-2)-3*H(j-3)+(K*(25*F(j+1)-48*F(j)+36*F(j-1)-16*F(j-2)+3*F(j-3))+12*hstep*(F(j)+1.5*tau+2*K)*h(j))/(1+0.5*K))/25;
    Theta(j+1)=((48*hstep-12*(hstep^2)*Pr*(f(j)-0.5*tau*hstep*(j+2))/(1+4*R/3))*Theta(j)-36*hstep*Theta(j-1)+16*hstep*Theta(j-2)-3*hstep*Theta(j-3)-12*Pr*Ec*(1+K)*((25*F(j+1)/12-4*F(j)+3*F(j-1)-4*F(j-2)/3+F(j-3)/4)^2)/(1+4*R/3)+12*(hstep^2)*(Pr*theta(j)*(F(j)+tau)-Pr*M*Ec*(F(j)^2)-A1*F(j)-B1*theta(j))/(1+4*R/3))/(25*hstep);
    Phi(j+1)=((48+12*hstep*tau*Sc*hstep*(j+2)-12*hstep*Sc*f(j))*Phi(j)-36*Phi(j-1)+16*Phi(j-2)-3*Phi(j-3)+12*hstep*phi(j)*Sc*(tau+F(j)))/25;
end
D=optimvar('D',4);   %% Lengendre Constants C
xi=exp(-1*eps);      %% From serth's paper the transformation from eta to xi
%syms f(xi);
fun=xi-1-log(xi)+power(log(xi),2).*(D(1)+D(2)*(2.*xi-1)+D(3)*(6*power(xi,2)-6.*xi+1)+D(4)*(20*power(xi,3)-30*power(xi,2)+12.*xi-1));
%fun=xi.*(1-(1./xi)+2*(xi-1).*(D(1)+D(2)*(2.*xi-1)+D(3)*(6*power(xi,2)-6.*xi+1)+D(4)*(20*power(xi,3)-30*power(xi,2)+12.*xi-1))+power(xi-1,2).*(2*D(2)+D(3)*(12.*xi-6)+D(4)*(60*power(xi,2)-60.*xi+12)));
obj=sum((fun-theta).^2);    %% Objective function
lsqproblem=optimproblem("Objective",obj); 
x0.D=[3/2 1/2 3/2 1/2];
[sol,fval] = solve(lsqproblem,x0);
%err=evaluate(obj,sol);
err=abs(Theta(N));
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
l3=l3-double(((l(exp(-ninf))-1)/m(exp(-ninf))));
Fi=smooth(F);
figure(1);
plot(eps,Fi,'b');
figure(2)
plot(eps,h);
figure(3)
plot(eps,theta);
figure(4)
plot(eps,phi);