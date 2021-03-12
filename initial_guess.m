function main
    xa = 0; xb =10;
    l1=-0.8;
    l2=1;
    l3=1;
    l4=1;
    solinit = bvpinit(linspace(xa,xb,100),[1 0.9 l1 -0.5*l1 l2 1 l3 1 l4]);
    options = bvpset('stat','on','RelTol',1e-6);
    sol=bvp4c(@bvpfcn,@bcfcn,solinit,options);
    xint = linspace(xa,xb,100);
    Sxint = deval(sol,xint);
    %plot(sol.y);
    plot(xint,Sxint(2,:));
end
function res=bcfcn(ya,yb)
S=1;
lambda=0.9;
res=[ya(1)-S;
    ya(2)-lambda;
    ya(4)+0.5*ya(3);
    ya(6)-1;
    ya(8)-1;
    yb(2);
    yb(4);
    yb(6);
    yb(8)];
end

function dydx =bvpfcn(x,y)
Sc=0.3;
tau=0.2;
delta=0.2;
M=0.5;
R=0.5;
Ec=0.3;
K=0.5;
A1=0.2;
B1=0.2;
delta1=0.2;
Pr=3.0;
dydx=[y(2);
    y(3);
    -1/(1+K)*(y(1)*y(3)-y(1)^2-tau*(y(1)+0.5*x*y(3))+K*y(5)+delta*y(6)+delta1*y(8)-M*y(2));
    y(5);
    -(1+0.5*K)*(y(1)*y(5)-y(2)*y(4)-tau*(1.5*y(4)+0.5*x*y(5))-K*(2*y(4)+y(3)));
    y(7);
    -(1+4*R/3)*(Pr*(y(1)*y(7)-y(2)*y(6))-Pr*tau*(2*y(6)+0.5*x*y(7))+Pr*Ec*(1+K)*(y(3)^2)+Pr*M*Ec*(y(2)^2)+A1*y(2)+B1*y(6));
    y(9);
    -(Sc*(y(1)*y(9)-y(2)*y(8))-Sc*tau*(y(8)+0.5*x*y(9)))];
end
