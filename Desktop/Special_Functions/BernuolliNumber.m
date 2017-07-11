% Bernuolli Number Calculation

function const=BernuolliNumber(N)

B(1)=1;
for n=2:N+2
    add=0;
    for k=1:n-1;
        C(n,1)=1;
        G(n,k)=(n-(k-1))./(k);
        C(n,k+1)= G(n,k).*C(n,k);
        add= add + C(n,k).*B(k);
    end
    if n>=4 && ~mod(n,2);
        B(n)=0;
    else
        B(n)=-add./C(n,k+1);
    end
end
B=B(B~=0);
if N>2 && mod(N,2);
    error('N should be 0, 1 or a positive even number!');
elseif N>=2 && ~mod(N,2);
    const=B(end);
else
    const=B(N+1);
end