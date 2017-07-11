%% S.Aksimsek, 2011
% lngamma_z Function

function lngamma_z=lngamma_z(rez,imz)

%n=12
if rez<6;
    recursion=ceil(6-rez); % Recursion number
    numbr=recursion + rez; % The value which will be begun from this value.
    if (rez<=0 && ceil(rez)==floor(rez)) % if the number rez<=0 and an integer,
                                    % "if (rez<=0 || ceil(rez)==floor(rez))"
                                    % (rez<=0 ||  round(rez)==rez)
        lngamma_z=NaN;               % lngamma_z_asymptotic_formula will not be a number
    else
        %         lngammaz=lngamma_z_asymptotic_formula(n,numbr,imz);
        while recursion>0
            numbr=numbr-1;
            numbr2=complex(numbr,imz);
            lngamma_z=lngamma_z-log(numbr2);
            recursion=recursion-1;
        end
    end
else
    lngamma_z=lngamma_z_asymptotic_formula(n,rez,imz);
end 
lngamma_z