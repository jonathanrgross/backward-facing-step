function [u0,v0] = IC_0(nx,ny)


u0 = zeros(ny,nx);
v0 = zeros(ny,nx);
% initial guess for u
y=0:1/(ny-1):1;
y=y';
k=1;
j=1;
    a=1.5;                                 	% height
    b=(round(nx/2)-round(nx/2)*j/nx)/nx;   	% location
    c=(round(nx/2)+round(nx/2)*j/nx)/nx;   	% width
    u0LBC = 4*a*((y-b)/c-((y-b)/c).^2);
    
    for k = 1:ny
        if u0LBC(k)<0
            u0LBC(k)=0;
        end
    end
u0(:,1) = u0LBC;
end


