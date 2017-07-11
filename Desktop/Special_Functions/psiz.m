%% S.Aksimsek, 2011
% Psi_z Function

function psiz=psiz(rez,imz)

n=14
if rez<6;
    recursion=ceil(6-rez); % Recursion number
    numbr=recursion + rez; % The value which will be begun from this value.
    if (rez<=0 && ceil(rez)==floor(rez)) % if the number rez<=0 and an integer,psiz will not be a number.
        psiz=NaN
    else
        psiz=psi_z_asymptotic_formula(n,numbr,imz)
        while recursion>0
            numbr=numbr-1;
            numbr2=complex(numbr,imz);
            psiz=psiz-1/(numbr2);
            recursion=recursion-1
        end
    end
else
    psiz=psi_z_asymptotic_formula(n,rez,imz);
end


