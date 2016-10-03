function B=banfac(B,N)
%
% programer: Tai Merzel
% AME 535
%
% This subroutine factorizes the band matrix from linear
% elements and puts them into LU form.  It uses the same type of setup
% as BANFAC.for given in the book.

% The tridiagonal system is solved

NP = N - 1;
for C = 1:NP
    CP = C + 1;
    B(2,CP) = B(2,CP)/B(3,C);
    B(3,CP) = B(3,CP) - B(2,CP)*B(4,C);
end
end % function banfac
