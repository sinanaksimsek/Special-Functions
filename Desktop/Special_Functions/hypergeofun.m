function [output error]= hypergeofun(alfa,beta,gama,z)
Nmax=1000; %maximum default value for N
n=0;
m=1;
sum=0;
term(m)= 1;
error=0;
%Calculation of minimum n value
a=beta*alfa*z;
b= alfa*z + beta*z-gama;
c=z-1;
x1= -0.5*(b/a) + sqrt( (0.5*(b/a))^2 - (z-1)/a );
x2= -0.5*(b/a) - sqrt( (0.5*(b/a))^2 - (z-1)/a );
if x1>0 && x2>0
    Nmin= ceil(1/ min(x1,x2));
elseif x1>0 && x2<0
    Nmin=ceil(1/x1);
elseif x2>0 && x1<0
    Nmin=ceil(1/x2);
else
    Nmin=1;
end
for n=1:Nmax
    m=m+1;
    R=(alfa+n)*(beta+n)*z/((gama+n)*(n+1)); %the ratio of terms t(n+1)/t(n)
    term(m)=R*term(m-1);
    maxT=max(term); %maximum of terms
    if term(m)<eps*maxT && R<1 && n>Nmin
        n
        output=sum
        break
    end
    sum=sum+term(m);
    error=error + eps*term(m)
end