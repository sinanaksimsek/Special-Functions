%% S.Aksimsek, 2011
function nlimit_of_HypergeometricFunction=nlimit_of_HypergeometricFunction(alfa,beta,gama,z)

a=beta*alfa*z;
b= alfa*z + beta*z-gama;
c=z-1;
x1= -0.5*(b/a) + sqrt( (0.5*(b/a))^2 - (z-1)/a );
x2= -0.5*(b/a) - sqrt( (0.5*(b/a))^2 - (z-1)/a );
if x1>0 && x2>0
    nlimit_of_HypergeometricFunction= ceil(1/ min(x1,x2));
elseif x1>0 && x2<0
    nlimit_of_HypergeometricFunction=ceil(1/x1);
elseif x2>0 && x1<0
    nlimit_of_HypergeometricFunction=ceil(1/x2);
else
    nlimit_of_HypergeometricFunction=1;
end