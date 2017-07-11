%% S.Aksimsek, 2011
% z0 Value Calculation for Psi Function

for i=1:101
    B2k(i)=BernuolliNumber(2*i); % Bernoulli numbers
end
for n=5:100
    for z=1:100;
        addd=0;
        for k=1:n-1;
            D=B2k(k)/(2*k*(z^(2*k)));
            addd = addd + D;
        end
        psi_z = log(z)-1/(2*z) + addd; %psi_z calculation
        if abs((D/psi_z))<eps % Observe relative error
            z0(n)=z;
            break
        end
    end
end
figure
n=5:100
z0=z0(5:100)
plot(n,z0)
xlabel('z0 values for psi asymptotic formula')
ylabel('n values')
[(5:100)' z0(1:96)']