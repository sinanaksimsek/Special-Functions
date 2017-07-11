function BesselRecurrenceBYGoldstein(M,v,z)

%Recurrence Techniques by M. Goldstein and R. M. Thaler

% M=20;    % Uppper level of Bessel Function
% z=0.5;
M0=floor(max(M,abs(z)));
n=M0+2;
Y(M0)=0;
Y(M0+1)=1;
flag=0;
% By using direct recursion formula, determine N0 and N values.
while n>M0;                             % Increase v in each step.
    Y(n)=(2*(n-1)/z)*Y(n-1)-Y(n-2);    % Recursion formula
    if (abs(1/Y(n))<eps & flag==0);      % Determine N0
        N0=n
        flag=1;
    end
    if abs(1/Y(n))<eps^2 ;               % Determine N
        N=n
        break
    end
    n=n+1;
end

Y=Y(M0+1:N)'  % Smaller than M0+1, all Y~ values will be zero.

% By using inverse recursion formula, calculate Jtilde values.

% Lets calculate J tilda
Jtilde(N+1)=0;
Jtilde(N)=1;
n=N+1;
while n>2;
   Jtilde(n-2)=-Jtilde(n)+(2*(n-1)/z)*Jtilde(n-1);
    n=n-1;
end
Jtilde=Jtilde(N+1:-1:1);
Jtilde=Jtilde'
disp(Jtilde)

% 1=2^v*SUM(v+2m)*(J_(v+2m)/z)*z^(-v)(gamma(v+m)/m!)
% For v=0, this formula reduces to the formula given below:
% J~0 + SUM(m=1,N0){J~2m}=1
sum=0;
m=1;
if (v==0);
    while (m<N0+1);
        if 2*m<N+1;
            term=Jtilde(2*m);
        else
            term=0;
        end
        m=m+1;
        sum=sum+term;
    end
    A=Jtilde(1)+2*sum;
        
else     %(v is not zero)
    while (m<N0+1);
        if 2*m<N+1;
            term= (v+2*m)*Jtilde(2*m)*exp(lngammaz(v+m,0))/factorial(m);
        else
            term=0;
        end
        m=m+1;
        sum=sum+term;
    end
    A=(2^v)*z^(-v)*sum;
end
A

% We have A coefficient. At this point we can calculate Bessel Function
% using Jtilde and A coefficient.

for m=1:M+1;
    J(m)=Jtilde(m)/A; %// to test A coefficient with Bessel functions as J~=A*J;
end   
B=J';
