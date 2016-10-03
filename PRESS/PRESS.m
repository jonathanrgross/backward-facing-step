function [u2,v2,P,Pend] = PRESS(u1,v1,nx,ny,dx,dy,dt,A,jpvt,Pend)
% this function implicitly solves the poision equation that defines the 
% pressure field and uses the pressure to calculate velocity.


% preallocate RHS
RHS = zeros(ny,nx);

% construct a matrix of values for the right hand side
for j = 2:nx-1
    for k = 2:ny-1;
        RHS(k,j) = (u1(k,j)-u1(k,j-1))/(dt*dx) + (v1(k,j)-v1(k-1,j))/(dt*dy);
    end
end

% because P(ny,nx) is fixed, define adjacent points in terms of Pend
RHS(ny-1,nx) = RHS(ny-1,nx) - Pend/(dy^2);
RHS(ny,nx-1) = RHS(ny,nx-1) - Pend/(dx^2);
RHS = reshape(RHS,[nx*ny 1]);
n=nx*ny;
RHS = RHS(1:n-1);

% solve for the pressure field
P = zeros(ny,nx);
[P]=solve(n-1,A,jpvt,RHS);   

% turn the vector P into an array
P(n)=Pend;
i = 1;
for j = 1:nx
    for k = 1:ny
        P(k,j) = P(i);
        i= i + 1;
    end
end

% assign corner values to average adjacent points
P(1,1) = (P(1,2) + P(2,1))/2;
P(1,nx) = (P(1,nx-1) + P(2,nx))/2;
P(ny,1) = (P(ny,2) + P(ny-1,1))/2;
P(ny,nx) = (P(ny,nx-1) + P(ny-1,nx))/2;

% calculate Pend in preparation for next iteration
Pend = sum(P(2:ny-1,nx))/length(P(2:ny-1,nx));

% calculate u2 and v2
for j = 2:nx-1
    for k = 2:ny-1
        u2(k,j) = u1(k,j) -dt*(P(k,j+1)-P(k,j-1))/(2*dx);
        v2(k,j) = v1(k,j) -dt*(P(k+1,j)-P(k-1,j))/(2*dy);
    end
end

% assign boundary values to u2 and v2
    % bottom boundary
    u2(1,:) = 0;
    v2(1,:) = 0;
    % top boundary
    u2(ny,:) = 0;
    v2(ny,:) = 0;
    % left boundary
    u2(:,1) = BC_parabolic(nx,ny);
    v2(:,1) = 0;
    % right boundary
    u2(:,nx) = u2(:,nx-1);
    v2(:,nx) = 0;
end




