% Calculation of Bessel Function
function Bn=BesselFunction(v,z)

Bn=0;  % First value of sum
n=1;
k=0;
%%%%%%%%%%%%%%%%%%%%
% Lets calculate the first term of the sum
g=1;
if v==0;
    g=1;
else
    for V=1:v;
        g=g*V;
    end
end
term(1)=(z/2).^(v)/((1)*g);
Bn= Bn + term(1);
%%%%%%%%%%%%%%%%%

%%%%% Bessel Function Calculation %%%
if z==0
    Bn=term(1)
else
    while k>=0
        n=n+1;
        k=k+1;
        R=(-1)*((z/2)^2)/((k)*(v+k)); %the ratio of terms t(k+1)/t(k)
        term(n)=R*term(n-1);
        Bn=Bn+term(n) ;
        maxT=max(term);             %maximum value of terms
        if (abs(term(n))<abs(term(n-1)) && abs(term(n))<eps*abs(maxT)); % Termination Condition
            break
        end
    end
end


