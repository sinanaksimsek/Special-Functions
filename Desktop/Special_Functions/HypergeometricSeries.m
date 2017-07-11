%% S.Aksimsek, 2011
function HypergeometricSeries=HypergeometricSeries(a,b,g,z)

error1=0;
error2=0;

if abs(z)>0.5  % If abs(z) is bigger than 0.5, the following formula will work
    n=g-a-b-1
    %************** n is an integer!**************************************
    if (n<0 & round(n)==n & imag(n)==0) % First of all let's calculate the phi function.
        n=-n;
        sum1=0;
        sum2=0;
        k=0;
        error=0;
        
        %******************** PHI FUNCTION*********************************
        
        % The first sum of the phi function
        if n~=1;
            term(1)=1/((n-1)*(a-1)*(b-1)*z);   % First term at k=1
            sum1=term(1);
            for k=2:n-1;
                R=(-1)*k/((n-k)*(a-k)*(b-k)*z);
                term(k)=R*term(k-1);
                sum1=sum1+term(k);
                error=error+abs(term(k))*eps;
            end
        else
            sum1=0;
            error=0;
        end
        
        % The second sum of the phi function
        Nmax=nlimit_of_HypergeometricFunction(a,b,n,z)-1;
        m=1;
        k=0;
        t1(1)=1;   % The first term of the FIRST part of  the second sum
        t2(1)= log(z)+psi(a)-psi(a-n+1)+psi(b)-psi(b-n+1)-psi(1); % The first term of the SECOND part of  the second sum
        term(1)=t1(1)*t2(1); % First term of the first sum at k=0
        sum2=term(1);
        for k=1:Nmax;
            m=m+1;
            R1=(a+k-1)*(b+k-1)*z/((k)*(n+k-1));
            t1(m)=R1*t1(m-1);
            t2(m)=log(z)+psi(a+k)-psi(a-n+1)+psi(b+k)-psi(b-n+1)+psi(n)-psi(n+k)-psi(k+1);
            term(m)=t1(m)*t2(m)
            sum2=sum2+term(m);
            error=error+abs(t1(m))*eps;   % BURDAN GELEN HATA term(m)'den GELÝYOR OLABÝLÝR!!!!!
        end        
        phi=sum1+sum2   % Phi function would be this.
        
        % *****************************************************************
        % **************** F(a,b,a+b-n+1,1-z) Equation ********************
        
        %Lets calculate F(a,b,n,z)
        m=1
        error1=0
        term(1)=1 % at k=0
        sum_f=term(1);
        for k=1:Nmax;
            m=m+1;
            R=(a+k-1)*(b+k-1)*z/((k)*(n+k-1));
            term(m)=R*term(m-1);
            sum_f=sum_f+term(m);
            error1 = error1 + abs(term(m))*eps; 
        end
        
        Part1=(((-1)^n)*exp(lngammaz(a+b-n+1,0)))/(exp(lngammaz(a-n+1,0)*exp(lngammaz(b-n+1,0))*(factorial(n-1))))
        Part2=(psi(a-n+1)+psi(b-n+1)-psi(n))*sum_f +phi
        error1=error1*Part2 +error
        FF=Part1*Part2;
        HypergeometricSeries=FF
        
    else  % n is not an integer
        Nmax=nlimit_of_HypergeometricFunction(a,b,a+b-g+1,1-z)-1;
        m=1;
        term1(1)=1;
        Sum1=term1(1);  % First sum at k=0
        for k=1:Nmax
            m=m+1;
            R1=(a+k-1)*(b+k)*(1-z)/((a+b-g+k)*(k));
            term1(m)=R1*term1(m-1);
            Sum1=Sum1+term1(m);
            error1=error1+abs(term1(m))*eps;
        end
        
        Nmax=nlimit_of_HypergeometricFunction(g-a,g-b,g-a-b+1,1-z)-1;
        m=1;
        term2(1)=1;
        Sum2=term2(1);  % First sum at k=0
        for k=1:Nmax;
            m=m+1;
            R2=(g-a+k-1)*(g-b+k-1)*(1-z)/((g-a-b+k)*(k));
            term2(m)=R2*term2(m-1);
            Sum2=Sum2+term2(m);
            error2=error2+abs(term2(m))*eps;
        end
        
        g1=exp(lngammaz(g,0));
        g2=exp(lngammaz(g-a-b,0));
        g3=exp(lngammaz(g-a,0));
        g4=exp(lngammaz(g-b,0));
        g5=exp(lngammaz(a+b-g,0));
        g6=exp(lngammaz(a,0));
        g7=exp(lngammaz(b,0));
        z_order=(1-z)^(g-a-b);
        
        error1=error1*(g1*g2/(g3*g4))
        error2=error2*(g1*g5/(g6*g7))*z_order
        error= error1+error2
        
        F=(g1*g2/(g3*g4))*Sum1 + (g1*g5/(g6*g7))*z_order*Sum2;
        HypergeometricSeries=F;
    end
    
else  % ( if abs(z)<0.5)
     Nmax=nlimit_of_HypergeometricFunction(a,b,g,z)-1;
     error=0
     m=1
     term(1)=1;
     Sum=term(1);  % First sum at k=0
     for k=1:Nmax;
            m=m+1
            R=(a+k-1)*(b+k-1)*z/((k)*(g+k-1));
            term(m)=R*term(m-1);
            Sum=Sum+term(m);
            error=error+abs(term(m))*eps;
     end
     HypergeometricSeries=Sum;
end

    