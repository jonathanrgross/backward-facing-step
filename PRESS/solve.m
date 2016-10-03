function [B]=solve(N,A,JPVT,B)
% author Tai Merzel
% AME 535
% Last revision: 12/1/2005
% SOLVE.m

% This program is simply a matlab version of SOLVE.for written
% by Fletcher; no major changes have been made.


NM = N - 1;
for K = 1:NM
    KP = K + 1;
    L = JPVT(K);
    S = B(L);
    B(L) = B(K);
    B(K) = S;
    for I = KP:N
        B(I) = B(I) + A(I,K)*S;
    end
end

for KA = 1:NM
    KM = N - KA;
    K = KM + 1;
    B(K) = B(K)/A(K,K);
    S = - B(K);
    for I = 1:KM
        B(I) = B(I) + A(I,K)*S;
    end
end
B(1) = B(1)/A(1,1);
end % function solve
