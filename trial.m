ninf=10;            %% n_infinity
h=0.02;             %% step size
n=(ninf)/h+3;      %% number of iterations
eta=-2*h:h:ninf;   %% Mesh Points
%xi=exp(-1*eta);      %% From serth's paper the transformation from eta to xi
syms xi;
D=[2.193002960090311, 4.784408827661294, 3.968568372188999, 1.441342999728768];
fun=xi-1-log(xi)+power(xi-1,2).*(D(1)+D(2)*(2.*xi-1)+D(3)*(6*power(xi,2)-6.*xi+1)+D(4)*(20*power(xi,3)-30*power(xi,2)+12.*xi-1));
%s=1+2*(D(1)+D(2)+D(3)+D(4));
g(xi)=diff(diff(fun));
l(xi)=diff(diff(diff(fun)));
s=1.2;
s=s+double((g(1)/l(1)))
x=[0.2 2];
fplot(g,x);