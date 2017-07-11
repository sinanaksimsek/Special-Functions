% S.Aksimsek, 2011 
%Miller Type Algorithm for
  
clear all
z=0.5;
V=20;
J(1)=BesselFunction(0,z)   %J0(z)
J(2)=BesselFunction(1,z)    %J1(z)
Bj(1)=besselj(0,z)
Bj(2)=besselj(1,z)
ErrorJ(1)=J(1)-Bj(1)
ErrorJ(2)=J(2)-Bj(2)
m=2
for v=1:V
    v
    J(m+1)=(2*(v)/z)*J(m)-J(m-1);
    Bj(m+1)=besselj(v+1,z);
    ErrorJ=Bj-J;
    m=m+1
end
T1=[J' Bj' ErrorJ']


%%%Apply Neumann Function to the recursion formula given below.


Y(1)=NeumannFunction(0,z);    %Y0(z)
Y(2)=NeumannFunction(1,z);    %Y1(z)
By(1)=bessely(0,z);
By(2)=bessely(1,z);
ErrorY(1)=Y(1)-By(1);
ErrorY(2)=Y(2)-By(2);

m=2
for v=1:V;
    Y(m+1)=(2*(v)/z)*Y(m)-Y(m-1);
    By(m+1)=bessely(v+1,z);
    ErrorY=By-Y;
    m=m+1;
end
T2=[Y' By' ErrorY']

