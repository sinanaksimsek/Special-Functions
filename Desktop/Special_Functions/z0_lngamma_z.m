%% S.Aksimsek, 2011
% z0 Value Calculation for the Asymptotic Formula of ln(gamma_z)

for i=1:101
    B2k(i)=BernuolliNumber(2*i); % Bernoulli numbers
end
for n=5:100
    for z=1:100;
        add=0;
        for k=1:n-1;
            D=B2k(k)/(2*k*(2*k-1)*(z^(2*k-1)));
            add = add + D;
        end
        lngamma_z = (z-1/2)*log(z)-z + 1/2*log(2*pi) + add; %lngamma_z calculation
        if abs((D/lngamma_z))<eps % Observe relative error
            z0(n)=z;
            break
        end
    end
end
figure
n=5:100
z0=z0(5:100)
plot(n,z0)
xlabel('z0 values for lngamma_z asymptotic formula')
ylabel('n values')
[(5:100)' z0(1:96)'] % List the z0 values for each n value on commond window.