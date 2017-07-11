%% S.Aksimsek, 2011
%%% Error representation for Bessel Function %%%
v=0; 
%v=20
m=1
for z=1:10:101
    B1 = BesselFunction(v,z); 
    B2 = besselj(v,z); 		 
    B3=B1-B2;
    BB(m,:)=[B1 B2 B3];
    m=m+1;
end
