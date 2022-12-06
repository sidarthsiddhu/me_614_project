! This program mimics staggered_grid_ins.m and employs PETSc to solve the linear system

Program Main

    implicit none
	
	integer :: nx,ny,nstep,maxIt,plot_freq
	Real :: Lx,Ly,Re,u_north,u_south,v_east,v_west
	Real :: time,dx,dy,dt,maxError,omega,tmp2,error
	
	integer :: is,i,j,it
	
	Real,Allocatable,Dimension(:,:) :: u,v,p,ut,vt,pold,uplot,vplot,tmp1
	Real,Allocatable,Dimension(:) :: x,y
	
! Set-up the problem constants and parameters
	Lx = 1.0
	Ly = 1.0
	Re = 100.0
	
	u_north = 1.0
	u_south = 0.0
	v_east = 0.0
	v_west = 0.0
	
	time = 0.0
	plot_freq = 100
	
	nx = 500
	ny = 500
	dx = Lx/nx
	dy = Ly/ny
	
	dt = 6e-04
	nstep = 50
	
	maxIt = 500
	maxError = 0.1
	omega = 1.85
	
! Allocate the main and supporting variables
	Allocate(u(nx+1,ny+2))
	Allocate(ut(nx+1,ny+2))
	Allocate(uplot(nx+1,ny+1))
	
	Allocate(v(nx+2,ny+1))
	Allocate(vt(nx+2,ny+1))
	Allocate(vplot(nx+1,ny+1))
	
	Allocate(p(nx+2,ny+2))
	Allocate(pold(nx+2,ny+2))
	Allocate(tmp1(nx+2,ny+2))
	
	Allocate(x(nx+1))
	Allocate(y(ny+1))
	
! Set-up x and y co-oridnates
	Do i=1,nx+1
		x(i) = (i-1)*(dx)
	Enddo
	
	Do j=1,ny+1
		y(j) = (j-1)*dy
	Enddo
	
! Set-initial values of u,v,p and support elements
	u = 0.0
	ut = 0.0
	 
	v = 0.0
	vt = 0.0
	 
	p = 0.0
	tmp1 = 0.0
	 
! Starting the time loop
	Do is=1,nstep
	
	! Set Tangential Velocity BC's
		u(:,1) = 2*u_south - u(:,2)
		u(:,ny+2) = 2*u_north - u(:,ny+1)
		v(1,:) = 2*v_west - v(2,:)
		v(nx+2,:) = 2*v_east - v(nx+1,:)
	
	! Calculate temproary u-velocity
		Do i=2,nx
			Do j=2,ny+1
				ut(i,j) = u(i,j) - dt*0.25/dx*( (u(i+1,j)+u(i,j))**2 -(u(i,j)+u(i-1,j))**2 )&
								&- dt*0.25/dx*( (u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j)) - (u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)) )&
								&+ dt/Re*( (u(i+1,j)-2*u(i,j)+u(i-1,j))/dx**2+ (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy**2)
			Enddo
		Enddo
		
		! Calculate temproary v-velocity
		Do i=2,nx+1
			Do j=2,ny
				vt = v(i,j) - dt*0.25/dx*( (u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j)) - (u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)) )&
						   &- dt*0.25/dy*( (v(i,j+1)+v(i,j))**2-(v(i,j)+v(i,j-1))**2 )&
						   &+ dt/Re*( (v(i+1,j)-2*v(i,j)+v(i-1,j))/dx**2+(v(i,j+1)-2*v(i,j)+v(i,j-1))/dy**2 )
		    Enddo
		Enddo
		
		Do i=2,nx+1
			Do j=2,ny+1
				tmp1(i,j) = 0.5/dt*( (ut(i,j)-ut(i-1,j))/dx + (vt(i,j)-vt(i,j-1))/dy )
			Enddo
		Enddo
		
		tmp2 = 1/(1/dx/dx + 1/dy/dy)
		
	! Solve for pressure using SOR
		Do it=1,maxIt
		
		! Set the old pressure to new
			pold = p
		
		! Set pressure ghost points 
			p(1,:) = p(2,:)
			p(nx+2,:) = p(nx+1,:)
			p(:,1) = p(:,2)
			p(:,ny+2) = p(:,ny+1)
			
			Do i=2,nx+1
				Do j=2,ny+1
					p(i,j) = (1.0-omega)*p(i,j) + omega*tmp2*( (.5/dx**2)*( p(i+1,j) + p(i-1,j) )+ (.5/dy**2)*( p(i,j+1) + p(i,j-1) ) - tmp1(i,j) )
				Enddo
			Enddo
		
		! Calculate error to stop iterations
			error = 0.0
			Do i=1,nx+2
				Do j=1,ny+2
					if(abs(p(i,j) - pold(i,j)) > error) error = abs(p(i,j) - pold(i,j))
				Enddo
			Enddo
			if( error < maxError) EXIT
			
		Enddo
		
	! Correct the u-velocity
		Do i=2,nx
			Do j=2,ny+1
				u(i,j) = ut(i,j) - dt/dx *(p(i+1,j) - p(i,j))
			Enddo
		Enddo
		
	! Correct the v-velocity
		Do i=2,nx+1
			Do j=2,ny
				v(i,j) = vt(i,j) - dt/dy *(p(i,j+1) - p(i,j))
			Enddo
		Enddo
		
		time = time + dt
		
	Enddo

! De-Allocate the main and supporitng variables
	Deallocate(u)
	Deallocate(ut)
	Deallocate(uplot)
	
	Deallocate(v)
	Deallocate(vt)
	Deallocate(vplot)
	
	Deallocate(p)
	Deallocate(pold)
	Deallocate(tmp1)
	
	Deallocate(x)
	Deallocate(y)

End Program
