%% S.Aksimsek, 2011
% Calculation of Hankel Function

function Hn=HankelFunction(v,kind,z)

if ceil(v)==floor(v);%If v is integer
    if kind==1;
        sgn=1;
    elseif kind==2
        sgn=-1;
    else
        'kind must be equal 1 or 2'
        return
    end
    Hn= BesselFunction(v,z) + sgn*i*NeumannFunction(v,z);
else
    %v is not an integer
    if kind==1;
        Hn= (BesselFunction(-v,z)- exp(-i*pi*v)*BesselFunction(v,z) ) / (i*sin(pi*v)) ;
    elseif kind==2;
        Hn= (exp(i*pi*v)*BesselFunction(v,z)-BesselFunction(-v,z) ) / (i*sin(pi*v)) ;
    else
        'Kind must be equal 1 or 2'
        return
    end
end
