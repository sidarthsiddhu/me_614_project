clear all; 
close all;

debug_version = 0;
print_result = 0;

Lx = 1.0; Ly = 1.0; Re=100;               % domain size and reynolds number
unorth=1.0; usouth=0; veast=0; vwest=0;   % boundary conditions

time=0.0; plot_freq = 100;

nx=250; ny=250; dx=Lx/nx; dy=Ly/ny; dt=0.0006;            % discretization   
nstep=6000; maxit=500; maxError=0.1; omg=1.85;

u  = zeros(nx+1,ny+2); ut   = u; uplot = zeros(nx+1,ny+1); % x-vel arrays
v  = zeros(nx+2,ny+1); vt   = v; vplot = zeros(nx+1,ny+1); % y-vel arrays
p  = zeros(nx+2,ny+2); tmp1 = p;                       

x = linspace(0,Lx,nx+1); y = linspace(0,Ly,ny+1);

tStart = tic;

for is=1:nstep

 u(1:nx+1,1) = 2*usouth-u(1:nx+1,2); u(1:nx+1,ny+2) = 2*unorth-u(1:nx+1,ny+1); % tangential vel BC
 v(1,1:ny+1) = 2*vwest -v(2,1:ny+1); v(nx+2,1:ny+1) = 2*veast -v(nx+1,1:ny+1); % tangential vel BC

 for i=2:nx; for j=2:ny+1      % temporary u-velocity (boundary values are not touched)
  ut(i,j) = u(i,j) - dt*0.25/dx*( (u(i+1,j)+u(i,j))^2 -(u(i,j)+u(i-1,j))^2 ) ...
                   - dt*0.25/dy*( (u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j)) - (u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)) ) ...
                   + dt/Re*( (u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2+ (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2) ;
 end; end

 for i=2:nx+1; for j=2:ny       % temporary v-velocity (boundary values are not touched)
  vt(i,j) = v(i,j) - dt*0.25/dx*( (u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j)) - (u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)) )...                                 
                   - dt*0.25/dy*( (v(i,j+1)+v(i,j))^2-(v(i,j)+v(i,j-1))^2 )...
                   + dt/Re*( (v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2+(v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2 ) ;    
 end; end   

% fprintf("----- Ut-Velocity --------\n")
% for i=1:nx+1
%     for j=1:ny+2
%         fprintf("%.16f ",ut(i,j));
%     end
%     fprintf("\n");
% end
% fprintf("\n");
% 
% fprintf("----- Vt-Velocity --------\n")
% for i=1:nx+2
%     for j=1:ny+1
%         fprintf("%.16f ",vt(i,j));
%     end
%     fprintf("\n");
% end
% fprintf("\n");

 for i=2:nx+1; for j=2:ny+1
  tmp1(i,j) = 0.5/dt*( (ut(i,j)-ut(i-1,j))/dx + (vt(i,j)-vt(i,j-1))/dy );
 end; end

%  fprintf("----- TMP1 --------\n")
% for i=1:nx+2
%     for j=1:ny+2
%         fprintf("%.16f ",tmp1(i,j));
%     end
%     fprintf("\n");
% end
% fprintf("\n");

 tmp2 = 1/(1/dx/dx + 1/dy/dy);

 for it=1:maxit	            % solve for pressure by SOR
   pold   = p;
   p(1,:) = p(2,:); p(nx+2,:) = p(nx+1,:); p(:,1) = p(:,2); p(:,ny+2) = p(:,ny+1); % set gosht values
   for i=2:nx+1; for j=2:ny+1
    p(i,j) = (1.0-omg)*p(i,j) ...
           + omg*tmp2*( (.5/dx^2)*( p(i+1,j) + p(i-1,j) )+ (.5/dy^2)*( p(i,j+1) + p(i,j-1) ) - tmp1(i,j) );
   end; end
   if max(max(abs(pold-p))) < maxError, break, end
 end
                                      
 for i=2:nx; for j=2:ny+1   % correct the u-velocity 
  u(i,j) = ut(i,j) - dt/dx*( p(i+1,j) - p(i,j) );
 end; end
      
 for i=2:nx+1; for j=2:ny   % correct the v-velocity
  v(i,j) = vt(i,j) - dt/dy*( p(i,j+1) - p(i,j) );
 end; end

if(debug_version == 1)
    fprintf("The number of iteration for nstep %d is: %d\n",is,it);
    if((is == nstep) && (print_result == 1))
        fprintf("----- Pressure --------\n")
        for i=1:nx+2
            for j=1:ny+2
                fprintf("%.16f ",p(i,j));
            end
            fprintf("\n");
        end
        fprintf("\n");
    end
end

% fprintf("----- U-Velocity --------\n")
% for i=1:nx+1
%     for j=1:ny+2
%         fprintf("%.16f ",u(i,j));
%     end
%     fprintf("\n");
% end
% fprintf("\n");
% 
% fprintf("----- V-Velocity --------\n")
% for i=1:nx+2
%     for j=1:ny+1
%         fprintf("%.16f ",v(i,j));
%     end
%     fprintf("\n");
% end
% fprintf("\n");

 time=time+dt;              
   
end % End of time step

%End time calculation
tEnd = toc(tStart);

%Plot final result
uplot(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1)); % uplot stores x-velocity at all points element corners 
vplot(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1)); % vplot stores y-velocity at all points element corners
figure(1); 
vm = sqrt(uplot'.*uplot' + vplot'.*vplot'); 
surf(x,y,vm,'facecolor','interp','edgecolor','none','facelighting','phong'); title('final velocity contour (MATLAB)'); colormap jet; colorbar; view([0,0,1]);axis tight;drawnow

fprintf("The time taken is: %.16f\n",tEnd);
fprintf("\n\n");

