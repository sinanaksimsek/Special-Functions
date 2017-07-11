% Binomial Coefficients Calculation

function BinomialCoefficient(N)
n=N;
C(n,1)=1;
for k=1:n;
if n>=2
G(n,k)=(n-(k-1))./(k);
C(n,k+1)=G(n,k).*C(n,k);
else
break; disp('n should be bigger than 1')
end
end
C(n,:)'