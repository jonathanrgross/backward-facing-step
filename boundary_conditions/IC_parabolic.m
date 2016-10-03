function [u0,v0] = IC_parabolic(nx,ny)


% initial guess for u
y=0:1/(ny-1):1;
y=y';
k=1;
for j = 1:nx
    a=1.5*(1-0.5*1/nx);                       % height
    b=(round(nx/2)-round(nx/2)*1/nx)/nx;      % location
    c=(round(nx/2)+round(nx/2)*1/nx)/nx;      % width
    u0(:,j) = 4*a*((y-b)/c-((y-b)/c).^2);
    
    for k = 1:ny
        if u0(k,j)<0
            u0(k,j)=0;
        end
    end
end

% initial guess for v
v0=zeros(ny,nx);
end
