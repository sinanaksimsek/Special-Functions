% Calculation of Chebyshev Polynomials by theta argument

function [T_theta U_theta]= ChebyshevPolynomial_theta(n,z)

teta=acos(z);
T_theta= cos(n*teta);
U_theta= sin((n+1)*teta)/sin(teta);
[T_theta U_theta]
