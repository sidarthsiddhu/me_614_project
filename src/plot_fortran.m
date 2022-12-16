%clear all;
close all;
Lx = 1.0;
Ly = 1.0;
uplotf = readmatrix("uplot.csv");
vplotf = readmatrix("vplot.csv");
pf = readmatrix("Pressure.csv");
vmf = sqrt(uplotf'.*uplotf' + vplotf'.*vplotf'); 
sz = size(vmf);
xf = linspace(0.0,Lx,sz(1));
yf = linspace(0.0,Ly,sz(2));
surf(xf,yf,vmf,'facecolor','interp','edgecolor','none','facelighting','phong'); title('final velocity contour (FORTRAN)'); colormap jet; colorbar; view([0,0,1]);axis tight;drawnow