%% S.Aksimsek, 2011
% lngamma_z Asymptotic Formula

function lngamma_z_asymptotic_formula=lngamma_z_asymptotic_formula(n, rez, imz);
z=complex(rez, imz);
for m=1:n
    B2k(m)=BernuolliNumber(2*m); % Bernoulli numbers
end
addd=0;
for k=1:n-1;
    D=B2k(k)/(2*k*(2*k-1)*(z^(2*k-1)));
    addd = addd + D;
end
lngamma_z_asymptotic_formula = (z-1/2)*log(z)-z + 1/2*log(2*pi) + addd; %lngamma_z calculation
