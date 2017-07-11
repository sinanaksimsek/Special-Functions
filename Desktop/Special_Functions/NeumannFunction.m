% Calculation of Neumann Function

function Yn=NeumannFunction(v,z)

if ceil(v)==floor(v);   % If v is an integer
    %************ The First Sum ************
    m=1;
    p=1;
    %**** (n-1)! calculation *****
    if v==0 && v==1;
        p=1;
    else
        for N=2:v;
            p=p*(N-1);
        end
    end
    %*****************************
    term1(1)=p*(z/2)^(-v);    % For k=0, the value of the first sum
    Sum1= term1(1);           % Firt value of Sum1
    if v==0;                    % If v=0, Sum1=0
        Sum1=0;
    elseif v==1;                % If v=1, Sum1=term1(1)
        Sum1=term1(1);
    else
        for k=1:v-1;           % If v>1, calculate Sum1
            m=m+1;
            R=1/((v-k)*(k))*((z/2)^2); %the ratio of terms t(k+1)/t(k)
            term1(m)=R*term1(m-1);
            Sum1=Sum1+term1(m);
        end
    end
    
    %*****************************************
    
    % ************ The Second Sum  *************
    k=0;
    g=1;
    m=1;
    
    %**** k! calculation *******************
    if v==0;
        g=1;
    else
        for V=1:v;
            g=g*V;
        end
    end
    %**************************************
    a(1)=((z/2)^v)/g;
    term2(1)= a(1)*(psi(v+1)+psi(1));   % For k=0, the value of the first sum
    Sum2=term2(1);                      % Firt value of Sum2
    if z==0;
        Sum2=term2(1);
        Yn= (1/pi)*(2*BesselFunction(v,z)*log(z/2)- Sum1-Sum2);
    else
        while k>=0;
            m=m+1;
            k=k+1;
            R=(-1)*(((z/2)^2)/((k)*(v+k)));  %the ratio of terms t(k+1)/t(k)
            a(m)=R*a(m-1);
            term2(m)=a(m)*(psi(v+k+1)+psi(k+1));
            Sum2=Sum2+term2(m);
            if (abs(term2(m))<abs(term2(m-1)) && abs(term2(m))<eps*abs(Sum2)); % Termination Condition
                break
            end
            
            %*****************************************
            Yn= (1/pi)*(2*BesselFunction(v,z)*log(z/2)- Sum1-Sum2);
        end
    end
else
    %v is not an integer
    Yn= (cos(pi*v)*BesselFunction(v,z)- BesselFunction(-v,z) ) /sin(pi*v);
end
