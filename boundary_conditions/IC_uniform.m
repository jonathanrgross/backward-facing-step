function [u0,v0] = IC_uniform(nx,ny)


u0 = 1.5*ones(ny,nx);
v0 = zeros(ny,nx);
u0(1:round(ny/2),1:nx)=0;

% boundary conditions:
    % u and v are 0 on the lower wall
    u0(1,:) = zeros(1,nx);
    v0(1,:) = zeros(1,nx);
    
    % u and v are 0 on the upper wall
    u0(ny,:) = zeros(1,nx);
    v0(ny,:) = zeros(1,nx);
    
    % u is zero on the lower half of left boundary
    u0(1:round(ny/2),1) = zeros(round(ny/2),1);
    
    % v is zero on the entire left boundary
    v0(1:round(ny/2),1) = zeros(round(ny/2),1);
end