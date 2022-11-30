! This program mimics staggered_grid_ins.m and employs PETSc to solve the linear system

Program Main

    implicit none
	
	integer :: nx,ny,nstep,maxIt,plot_freq
	Real :: Lx,Ly,rRe,u_north,u_south,v_east,v_west
	Real :: time,dx,dy,dt,maxError,omega,tmp2
	
	integer :: is,i,j,it
	
	Real,Allocatable,Dimension(:,:) :: u,v,p,ut,vt,p_old,uplot,vplot,tmp1
	Real,Allocatable,Dimension(:) :: x,y
	
	!Set-up the problem constants and parameters
	Lx = 1.0
	Ly = 1.0
	Re = 100.0
	
	u_north = 1.0
	u_south = 0.0
	v_east = 0.0
	v_west = 0.0
	
	time = 0.0
	plot_freq = 100
	
	nx = 100
	ny = 100
	dx = Lx/nx
	dy = Ly/ny
	
	dt = 6e-04
	nstep = 6000
	
	maxIt = 500
	maxError = 0.1
	omega = 1.85
	
	!Allocate the main and supporting variables
	Allocate(u(nx+1,ny+2))
	Allocate(ut(nx+1,ny+2))
	Allocate(uplot(nx+1,ny+1))
	
	Allocate(v(nx+2,ny+1))
	Allocate(vt(nx+2,ny+1))
	Allocate(vplot(nx+1,ny+1))
	
	Allocate(p(nx+2,ny+2))
	Allocate(p_old(nx+2,ny+2))
	Allocate(tmp1(nx+2,ny+2))
	
	Allocate(x(nx+1))
	Allocate(y(ny+1))
	
	!Set-up x and y co-oridnates
	Do i=1,nx+1
		x(i) = (i-1)*(dx)
	Enddo
	
	Do i=1,nx+2
		y(j) = (j-1)*dy
	Enddo
	
	!Set-initial values of u,v,p and support elements
	 u = 0.0
	 ut = 0.0
	 
	 v = 0.0
	 vt = 0.0
	 
	 p = 0.0
	 tmp1 = 0.0

	!De-Allocate the main and supporitng variables
	Deallocate(u)
	Deallocate(ut)
	Deallocate(uplot)
	
	Deallocate(v)
	Deallocate(vt)
	Deallocate(vplot)
	
	Deallocate(p)
	Deallocate(p_old)
	Deallocate(tmp1)
	
	Deallocate(x)
	Deallocate(y)

End Program
