N=25;
f=zeros(1,5);
h=zeros(1,5);
theta=zeros(1,5);
phi=zeros(1,5);
F=zeros(1,5);
H=zeros(1,5);
Theta=zeros(1,5);
Phi=zeros(1,5);
hstep=10/N;
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
figure(1)
plot(f);
figure(2)
plot(h);
figure(3)
plot(F);
figure(4)
plot(H);
%f(1)=1;
%f(2)=0.9;
%f(3)=-2.517;
%f(4)=1.258;
%f(5)=-2.228;

%i=5;
%F(i)=(f(i+1)/4+(5*f(i))/6-(3*f(i-1))/2+f(i-2)/2-f(i-3)/12)/hstep;
%H(i)=(h(i+1)/4+(5*h(i))/6-(3*h(i-1))/2+h(i-2)/2-h(i-3)/12)/hstep;
%Theta(i)=(theta(i+1)/4+(5*theta(i))/6-(3*theta(i-1))/2+theta(i-2)/2-theta(i-3)/12)/hstep;
%Phi(i)=(phi(i+1)/4+(5*phi(i))/6-(3*phi(i-1))/2+phi(i-2)/2-phi(i-3)/12)/hstep;