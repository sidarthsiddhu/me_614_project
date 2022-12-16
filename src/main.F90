! This program mimics staggered_grid_ins.m and employs PETSc to solve the linear system

Program Main

    implicit none
	
	integer :: nx,ny,nstep,maxIt,plot_freq
	Real :: Lx,Ly,Re,u_north,u_south,v_east,v_west
	Real :: time,dx,dy,dt,tol,omega,tmp2,errorMax,t1,t2
	
	integer :: is,i,j,it,debug_version,print_result
	
	Real(kind=8),Allocatable,Dimension(:,:) :: u,v,p,ut,vt,pold,uplot,vplot,tmp1,error
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
	
	nx = 200
	ny = 200
	dx = Lx/nx
	dy = Ly/ny
	
	dt = 0.0006
	nstep = 6000
	
	maxIt = 500
	tol = 0.01
	omega = 1.85
	
	debug_version = 0
	print_result = 0
	
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
	Allocate(error(nx+2,ny+2))
	
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
	
! start time
	Call CPU_TIME(t1)
	
! Starting the time loop
	Do is=1,nstep
	
	! Set Tangential Velocity BC's
		Do i=1,nx+1
			u(i,1) = 2*u_south - u(i,2)
			u(i,ny+2) = 2*u_north - u(i,ny+1)
		Enddo
		
		Do j=1,ny+1
			v(1,j) = 2*v_west - v(2,j)
			v(nx+2,j) = 2*v_east - v(nx+1,j)
		Enddo
		
	! Calculate temproary u-velocity
		Do i=2,nx
			Do j=2,ny+1
				ut(i,j) = u(i,j) - dt*(0.25/dx)*( (u(i+1,j)+u(i,j))**2 -(u(i,j)+u(i-1,j))**2 )&
								&- dt*(0.25/dx)*( (u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j)) - (u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)) )&
								&+ (dt/Re)*( (u(i+1,j)-2*u(i,j)+u(i-1,j))/(dx**2) + (u(i,j+1)-2*u(i,j)+u(i,j-1))/(dy**2))
			Enddo
		Enddo
		
		! Calculate temproary v-velocity
		Do i=2,nx+1
			Do j=2,ny
				vt(i,j) = v(i,j) - dt*0.25/dx*( (u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j)) - (u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)) )&
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
			error(:,:) = ABS(p(:,:) - pold(:,:))
			errorMax = MAXVAL(error)
			if(errorMax <= tol) EXIT
			
		Enddo
		
	! Correct the u-velocity
		Do i=2,nx
			Do j=2,ny+1
				u(i,j) = ut(i,j) - ((dt/dx) *(p(i+1,j) - p(i,j)))
			Enddo
		Enddo
		
	! Correct the v-velocity
		Do i=2,nx+1
			Do j=2,ny
				v(i,j) = vt(i,j) - ((dt/dy) *(p(i,j+1) - p(i,j)))
			Enddo
		Enddo
		
		time = time + dt
		
		if(debug_version .EQ. 1) then
			write(6,'(a,I0,a,I0)') "The number of iteration for nstep ",is," is: ",it
			if((is .EQ. nstep) .AND. (print_result .EQ. 1)) Call print_variable(p,nx+2,ny+2,1,nx+2,1,ny+2,"Pressure")
		Endif
		
	Enddo

! End time
	Call CPU_TIME(t2)
	
	Do i=1,nx+1
		Do j=1,ny+1
			uplot(i,j) = 0.5*(u(i,j+1) + u(i,j))
			vplot(i,j) = 0.5*(v(i+1,j) + v(i,j))
		Enddo
	Enddo

	write(6,'(a,F)') "The time taken to run the iterations are: ",t2-t1

! Write plot velocities out to plot it in matlab
	Call write_to_file(uplot,nx+1,ny+1,"uplot.csv")
	Call write_to_file(vplot,nx+1,ny+1,"vplot.csv")
	
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
	Deallocate(error)
	
	Deallocate(x)
	Deallocate(y)

End Program

Subroutine print_variable(k,npts_x,npts_y,x_start,x_end,y_start,y_end,k_name)

	Integer,intent(in) :: npts_x,npts_y,x_start,x_end,y_start,y_end
	Real(kind=8),dimension(npts_x,npts_y),intent(in) :: k
	Character(len=500),intent(in) :: k_name 
	Integer :: i,j
	
	write(6,'(a)') " "
	write(6,'(a,a11,a)') "---------- Printing  ",TRIM(k_name)," ----------"
	Do i=x_start,x_end
		Do j=y_start,y_end
			write(6,'(F,a)',advance="no") k(i,j)," "		
		Enddo
		write(6,'(a)') " "
	Enddo
	write(6,'(a,a11,a)') "---------- End Printing  ",TRIM(k_name)," ----------"
	write(6,'(a)') " "

End Subroutine

subroutine write_to_file(k,npts_x,npts_y,k_name)

	Integer,intent(in) :: npts_x,npts_y
	Real(kind=8),dimension(npts_x,npts_y),intent(in) :: k
	Character(len=500),intent(in) :: k_name
	Character(len=500) file_name
	
	Open(unit=8,file=TRIM(k_name),status='replace')
	Do i=1,npts_x
		Do j=1,npts_y
			If(j .NE. npts_y) then
				write(8,'(F,a)',advance="no") k(i,j),","
			Else
				write(8,'(F)') k(i,j)
			Endif
		Enddo
	Enddo
	Close(8)
	
End subroutine
	
	
	
	