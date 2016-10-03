function X=bansol(R,B,N)
% programmer Tai Merzel
% AME 535


% This subroutine uses LU factorization given R
% and solves for X.  It uses the same type of setup
% as BANSOL.for given in the book.

% The tridiagonal system is solved
A = N - 1;
for C = 1:A
    JP= C+1;
    R(JP)= R(JP)-B(2,JP)*R(C);
end  
X(N) = R(N)/B(3,N);
for C = 1:A;
    CA= N-C;
    X(CA)=(R(CA)-B(4,CA)*X(CA+1))/B(3,CA);
end
end % function bansol
