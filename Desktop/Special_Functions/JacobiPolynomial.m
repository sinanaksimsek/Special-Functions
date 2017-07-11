%% S.Aksimsek, 2011
%Calculation of Jacobi Polynomials by recursive formula

function output=JacobiPolynomial(N,z,alpha,beta)

% Initial values of P%%%%%%%%%%%%%%%%%%%%%%%%%%%
P(1)=1; %P0 at n=0
P(2)=0.5* (2*(alpha+1)+(alpha+beta+2)*(z-1)); %P1 at n=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m=2:N;
n=m-1;
alpha_n=(2*(n+1)*(n+alpha+beta+1))/((2*n+alpha+beta+1)*(2*n+alpha+beta+2));
beta_n=(beta^2-alpha^2)/((2*n+alpha+beta)*(2*n+alpha+beta+2));
gamma_n=(2*(n+alpha)*(n+beta))/((2*n+alpha+beta)*(2*n+alpha+beta+1));
P(m+1)=((z-beta_n)*P(m)-gamma_n*P(m-1))/alpha_n;
end
output=P(end);
