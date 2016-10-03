function [u4,v4]=VISC(u2,v2,Re,nx,ny,dx,dy,dt)
% this function solves for the viscous terms in the navier stokes equation


% preallocate
u3=u2;
v3=v2;

% first viscous step: diffusion in x direction
for k = 1:ny
    u3(k,:)=VISCneu(u2(k,:),nx,dx,dt,Re);
    v3(k,:)=VISCneu(v2(k,:),nx,dx,dt,Re);
end


% preallocate
u4=u3;
v4=v3;

% second viscous step: diffusion in y direction
for j = 1:nx
    u4(:,j)=VISCdir(u3(:,j),ny,dy,dt,Re);
    v4(:,j)=VISCdir(v3(:,j),ny,dy,dt,Re);
end

% reassign values to boundaries
u4(ny,:)=0;
u4(:,1) = BC_parabolic(nx,ny);
u4(:,nx) = u4(:,nx-1);
v4(ny,:)=0;
v4(:,nx)=0;
end
