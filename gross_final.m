function gross_final(Re,nx,ny,IC,tolerance)

%Re=40;
%nx = 61;
%ny = 41;
%IC = 1;
%tolerance = 0.01;


% declare variables
dx=6/(nx-1);
dy=2/(ny-1);
n=nx*ny;
dt=0.01;
tmax = 20;



% construct initial conditions
if IC == 1
    % uniform upper half flow
    [u0,v0] = IC_uniform(nx,ny);
elseif IC == 2
    % parabolic initial conditions
    [u0,v0] = IC_parabolic(nx,ny);
elseif IC == 3
    % parabolic initial conditions with interpolation
    [u0,v0] = IC_parabolic_iterpolate(nx,ny);
else
    % zero everywhere except at the boundary
    [u0,v0] = IC_0(nx,ny);
end


% calculate ideal outflow velocity profile
u_out = parabolic_outflow(nx,ny);


% construct the factorized matrix that will be used to solve for pressure
[A,jpvt]=constructA(nx,ny,dx,dy);


% initialize some variables
Fu_old = 0;
Fv_old = 0;
time = 0;
iteration = 1;
Pend = 1;
diff = ones(1,tmax/dt);



tic
while diff(iteration) > tolerance && time < tmax
    iteration = iteration + 1;
    
    % Nonlinear step
    [u1,v1,Fu_old,Fv_old]=CONVEC(u0,v0,Re,nx,ny,dx,dy,dt,Fu_old,Fv_old,iteration);
    
    % Pressure step
    [u2,v2,P,Pend] = PRESS(u1,v1,nx,ny,dx,dy,dt,A,jpvt,Pend);
    
    % Viscous step
    [u4,v4]=VISC(u2,v2,Re,nx,ny,dx,dy,dt);
    
    diff(iteration) = sum(sum(abs(u4-u0) + abs(v4-v0)));
    u_rms(iteration) = sqrt(sum((u_out-u4(:,nx)).^2/(ny-2)));
    Q_rms(iteration) = sqrt(sum((1 - 2*sum(u4)/ny).^2/nx));
    u0=u4;
    v0=v4;
    time = time + dt;
end
toc

ufinal=u4;
vfinal=v4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 4a. u profiles along x
figure
for j = 1:2:nx
    actual = plot(ufinal(:,j),-1:dy:1);
    hold on
end
u_out = parabolic_outflow(nx,ny);
ideal = plot(u_out,-1:dy:1,'r--');
xlabel('velocity')
ylabel('y')
title('u profiles')
legend([actual,ideal],'u profiles at various locations','ideal parabolic exit velocity','location','SouthEast')


% 4b. v profiles along x
figure
cmap=colormap(jet);
for j = 2:2:nx
    pick_color = cmap(round(60*j/nx),:);
    plot(vfinal(:,j),-1:dy:1,'Color',pick_color);
    hold on
end
xlabel('velocity')
ylabel('y')
title('v profiles')


% 5a. u along the centerline
figure
plot(0:dx:6,ufinal((ny+1)/2,:))
hold on
plot([0,6],[0.75, 0.75],'r--')
xlabel('x')
ylabel('u')
title('u vs x at centerline')
legend('u along centerline','centerline u of ideal outflow','location','SouthEast')


% 5b. v along the centerline
figure
plot(0:dx:6,vfinal((ny+1)/2,:))
hold on
plot([0,6],[0, 0],'r--')
xlabel('x')
ylabel('v')
title('v vs x at centerline')
legend('v along centerline','centerline v of ideal outflow','location','SouthEast')


% 6. flow rate
figure
Q = 2*sum(ufinal)/ny;
plot(0:dx:6,Q)
hold on
plot([0 6],[1 1],'r--')
xlabel('x')
ylabel('flowrate')
title('flow rate along the length of the domain')
axis([0 6 0 1.2*max(Q)])
legend('flowrate','flowrate of ideal outflow','location','SouthEast')



screensize = get( groot, 'Screensize');
% 1. velocity vector plot
figure('Position', [50, screensize(4)-400, 900, 300]);
[x,y] = meshgrid(0:dx:6,-1:dy:1);
quiver(x,y,ufinal,vfinal)
axis equal
title(['vector plot, Re = ', num2str(Re)])
hold on
plot([0 0 6],[0 -1 -1],'LineWidth',3,'Color','k')
plot([0 6],[1 1],'LineWidth',3,'Color','k')
hold on
[M,I]=min(abs(ufinal(2,10:nx)));
reattach_point = 6*(I+9)/nx;
plot(reattach_point,-1,'rx','MarkerSize',12,'LineWidth',2)


% 2. streamline plot
figure('Position', [500, 50, 900, 300]);
starty = -1:dy:1;
startx = starty*0;
streamline(x,y,ufinal,vfinal,startx,starty)
axis equal
title(['streamlines, Re = ', num2str(Re)])
hold on
plot([0 0 6],[0 -1 -1],'LineWidth',3,'Color','k')
plot([0 6],[1 1],'LineWidth',3,'Color','k')


% 3. reattachment length
hold on
[M,I]=min(abs(ufinal(2,10:nx)));
reattach_point = 6*(I+9)/nx;
plot(reattach_point,-1,'rx','MarkerSize',12,'LineWidth',2)


% output reattachment length and rms error values
disp(['reattachment length = ',num2str(reattach_point)])
disp(['u_rms = ',num2str(u_rms(end))])
disp(['Q_rms = ',num2str(Q_rms(end))])
disp(['iterations = ',num2str(iteration)])
end

