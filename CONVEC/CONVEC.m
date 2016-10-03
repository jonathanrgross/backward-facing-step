function [u1,v1,Fu_old,Fv_old]=CONVEC(u0,v0,Re,nx,ny,dx,dy,dt,Fu_old,Fv_old,iteration)
% This function explicitly solves for the nonlinear step using central
% differecing.  It solves only for interior points, and does not impart
% boundary conditions.


for j = 2:nx-1
    jm=j-1;
    jp=j+1;
    for k = 2:ny-1
        km=k-1;
        kp=k+1;

        Fu(k,j) = -u0(k,j)*(u0(k,jp)-u0(k,jm))/(2*dx) ...
            - v0(k,j)*(u0(kp,j)-u0(km,j))/(2*dy) + 3/(2*Re);
        
        Fv(k,j) = -u0(k,j)*(v0(k,jp)-v0(k,jm))/(2*dx) ...
            - v0(k,j)*(v0(kp,j)-v0(km,j))/(2*dy);
    end
end

% for the first iteration set F matricies to equal current iteration
if iteration-1 == 1
    Fu_old = Fu;
    Fv_old = Fv;
end


u1=u0;
v1=v0;

u1(2:ny-1,2:nx-1) = dt*((3/2)*Fu(2:ny-1,2:nx-1) - (1/2)*Fu_old(2:ny-1,2:nx-1)) + u0(2:ny-1,2:nx-1);
v1(2:ny-1,2:nx-1) = dt*((3/2)*Fv(2:ny-1,2:nx-1) - (1/2)*Fv_old(2:ny-1,2:nx-1)) + v0(2:ny-1,2:nx-1);

Fu_old = Fu;
Fv_old = Fv;
end


