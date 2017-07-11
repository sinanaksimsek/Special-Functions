% Calculation of Chebyshev Polynomials by z argument
 
function [T U]=ChebyshevPolynomial(n,z)
 
S1=exp(lngammaz(1/2+n,0))/ exp(lngammaz(1/2,0));
S2= exp(lngammaz(3/2+n,0))/ exp(lngammaz(3/2,0));
Jp1=JacobiPolynomial(n,z,-1/2,-1/2); %Jacobi Polynomial Calculation
Jp2=JacobiPolynomial(n,z,1/2,1/2);   %Jacobi Polynomial Calculation
T= (factorial(n)/S1)*Jp1;
U= (factorial(n+1)/S2)*Jp2;
[T U]

function [T U]=chebyfun(n,z)
% Calculation of Chebyshev Polynomials by theta argument
teta=acos(z);
T= cos(n*teta);
U= sin((n+1)*teta)/sin(teta);