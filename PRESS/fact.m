function [A,JPVT]=fact(N,A,JPVT)
% programer: Tai Merzel
% AME 535
% Last revision: 12/1/2005
% FACT.m

% This program is simply a matlab version of FACT.for written
% by Fletcher; no major changes have been made.

NM = N - 1;

for K = 1:NM
    KP = K + 1;
    
    L = K;
    for I = KP:N
        if(abs(A(I,K))>abs(A(L,K)))
            L = I;
        end
    end
    JPVT(K) = L;
    S = A(L,K);
    A(L,K) = A(K,K);
    A(K,K) = S;
    
    for I = KP:N
        A(I,K) = -A(I,K)/S;
    end
    
    for J = KP:N
        S = A(L,J);
        A(L,J) = A(K,J);
        A(K,J) = S;
        if (abs(S)>1.0E-15)  
            for I = KP:N
                A(I,J) = A(I,J) + A(I,K)*S;
            end
        end
    end
end
end % function fact
