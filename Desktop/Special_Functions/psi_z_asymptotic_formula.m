%% S.Aksimsek, 2011
% psi_z Asymptotic Formula

function psi_z_asymptotic_formula=psi_z_asymptotic_formula(n, rez, imz);
z=complex(rez, imz);
for m=1:n
    B2k(m)=BernuolliNumber(2*m); % Bernoulli numbers
end
addd=0;
for k=1:n-1;
    D=B2k(k)/(2*k*(z^(2*k)));
    addd = addd + D;
end
psi_z_asymptotic_formula = log(z)-1/(2*z) + addd; %psi_z calculation
