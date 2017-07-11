%% S.Aksimsek, 2011
%Calculation of Laguerre Polynomials by recursive formula

function output=LaguerrePolynomial(N,z)

alpha=0;

% Initial values of P%%%%%%%%%%%%%%%%%%%%%%%%%%%
L(1)=1; %PO at n=0
L(2)=-z+1; %P1 at n=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m=2:N
n=m-1;
alpha_n=-(n-1);
beta_n=2*n+alpha+1;
gamma_n=-(n+alpha);
L(m+1)=((z-beta_n)*L(m)-gamma_n*L(m-1))/alpha_n;
end
output=L(end);