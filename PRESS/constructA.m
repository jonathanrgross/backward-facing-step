function [A,jpvt]=constructA(nx,ny,dx,dy)
% This function generates the matrix of coefficients that will be be used
% to solve for pressure.  It first assigns values to five diagonals in a
% matrix, and then factorizes the matrix.

% define coefficients
n=nx*ny;
d1=1/dx^2;
d2=1/dy^2;
d3=-2/dx^2 - 2/dy^2;
d4=1/dy^2;
d5=1/dx^2;


%generate lowest diagonal
    ve5=d1*ones(1,n-ny);
    ve5(length(ve5)-ny+1:length(ve5))=2*ve5(length(ve5)-ny+1:length(ve5));
    D5 = diag(ve5,-ny);

    
%generate diagonal -1
    ve4=d2*ones(1,n-1);
    for i = 1:n+ny
        if mod(i,ny)==0
            if i<=n-1
            ve4(i)=0;
            end
            if i-1<=n-1
            ve4(i-1)=2*ve4(i-1);
            end
        end
    end
    D4 = diag(ve4,-1);

    
%generate main diagonal
    ve1=d3*ones(1,n);
    D1 = diag(ve1,0);

    
%generate diagonal +1
    ve2=d4*ones(1,n-1);
    for i = 1:n+ny
        if mod(i,ny)==0
            if i <=n-1
            ve2(i)=0;
            end
            if i <=n-1+ny
            ve2(i-ny+1)=2*ve2(i-ny+1);
            end
        end
    end
    D2 = diag(ve2,1);

    
%generate highest diagonal
    ve3=d5*ones(1,n-ny);
    ve3(1:ny)=2*ve3(1:ny);
    D3 = diag(ve3,ny);

    
% superimpose the diagonal matricies together to get the coefficient matrix
A = D1+D2+D3+D4+D5;

% remove last row and column, since P(ny,nx) is fixed
A = A(1:n-1,1:n-1);
jpvt = zeros(1,n-1);

% factorize the matrix using the function fact
[A,jpvt]=fact(n-1,A,jpvt);
end