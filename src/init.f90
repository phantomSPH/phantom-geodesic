module init
implicit none
real, parameter :: pi = acos(-1.)
contains

!--- Wrapper subroutine called in main.f90
subroutine setup(xall,vall,np)
   integer, intent(out) :: np
   real, allocatable, intent(inout), dimension(:,:) :: xall,vall
   integer, parameter :: isetup = 3

   select case(isetup)
   case(1)
      call setup_multisphere(xall,vall,np)
   case(2)
      call setup_2testbodies(xall,vall,np)
   case(3)
      call setup_singletype(xall,vall,np)
   case(4)
      call setup_multitypes(xall,vall,np)
   case(5)
      call setup_dude(xall,vall,np)
   case(6)
      call setup_sphere(xall,vall,np)
   end select

end subroutine setup

!--- Setup up a single body (wrapper to call initialise once)
subroutine setup_singletype(xall,vall,np)
   integer, intent(out) :: np
   real, allocatable, intent(inout), dimension(:,:) :: xall,vall
   np = 1
   allocate(xall(3,np),vall(3,np))
   call initialise(xall(:,1),vall(:,1),'precession')
end subroutine setup_singletype

!--- Setup up multipe (wrapper to call initialise multiple times)
subroutine setup_multitypes(xall,vall,np)
   integer, intent(out) :: np
   real, allocatable, intent(inout), dimension(:,:) :: xall,vall
   real :: r0
   integer :: i
   np = 10
   allocate(xall(3,np),vall(3,np))
   do i=1,np
      r0 = 0.+i*4.
      call initialise(xall(:,i),vall(:,i),'radial',r0)
   enddo
end subroutine setup_multitypes

!--- Setup a proxy for a sphere of particles
subroutine setup_sphere(xall,vall,np)
   integer, intent(out) :: np
   real, allocatable, intent(inout), dimension(:,:) :: xall,vall
   real :: dlayer, dr, rotate_y(3,3)
   integer :: i,j,k, nlayers, nrings,index,n,nringsmax
   real, parameter :: translate(3) = (/90.,0.,0./)
   dlayer = 0.5
   dr = 0.5
   nlayers = 7
   nringsmax = 4
   index = 0
   np=107
   allocate(xall(3,np),vall(3,np))
   do i=0,(nlayers-1)/2 !layer index
      nrings = nringsmax-i
      do j=1,nrings !number of rings in layer
         n=5*(j-1)+1
         do k=1,n !number of points in ring
            index=index+1
            if (j==0) then
               xall(:,index)= (/0.,0.,i*dlayer/)
            else
               xall(:,index)= (/(j-1)*dr*cos(2.*pi/n*k),(j-1)*dr*sin(2.*pi/n*k), i*dlayer/)

               index=index+1
               xall(:,index)= (/(j-1)*dr*cos(2.*pi/n*k),(j-1)*dr*sin(2.*pi/n*k),-i*dlayer/)

            endif
         enddo
      enddo
   enddo
   call get_rotation_matrix(-pi/4.,rotate_y,'y')
   do i=1,np
      xall(:,i)=xall(:,i) + translate
      xall(:,i)=matmul(rotate_y,xall(:,i))
      vall(:,i)=(/0.,0.0521157,0./)
   enddo
end subroutine setup_sphere

!--- Setup multiple spheres of particles
subroutine setup_multisphere(xall,vall,np)
   integer, intent(out) :: np
   real, allocatable, intent(inout), dimension(:,:) :: xall,vall
   real, allocatable, dimension(:,:) :: x1,v1,x2,v2,x3,v3
   integer :: n1, n2,i,n3
   real :: rotate_y(3,3), rotate_x(3,3), rotate_z(3,3),theta_y, theta_z, theta_x

   theta_y=-0.5*pi
   call get_rotation_matrix(theta_y,rotate_y,'y')

   theta_z=-0.5*pi
   call get_rotation_matrix(theta_z,rotate_z,'z')

   theta_x= 0.5*pi
   call get_rotation_matrix(theta_x,rotate_x,'x')

   call setup_sphere(x1,v1,n1)
   call setup_sphere(x2,v2,n2)
   call setup_sphere(x3,v3,n3)
   np = n1+n2+n3
   allocate(xall(3,np),vall(3,np))
   do i=1,n2
      x2(:,i)=matmul(rotate_y,x2(:,i))
      v2(:,i)=matmul(rotate_y,v2(:,i))
   enddo
   do i=1,n3
      x3(:,i)=matmul(rotate_x,matmul(rotate_z,x3(:,i)))
      v3(:,i)=matmul(rotate_x,matmul(rotate_z,v3(:,i)))
   enddo
   xall(:,1:n1) = x1
   vall(:,1:n1) = v1
   xall(:,n1+1:n1+n2) = x2
   vall(:,n1+1:n1+n2) = v2
   xall(:,n1+n2+1:n1+n2+n3) = x3
   vall(:,n1+n2+1:n1+n2+n3) = v3

end subroutine setup_multisphere

!--- Setup a stick-figure person
subroutine setup_dude(xall,vall,np)
   integer, intent(out) :: np
   real, allocatable, intent(inout), dimension(:,:) :: xall,vall
   integer :: i,index,nleg,nbody,narm,nhead
   real :: rotate_z(3,3),head_radius,rtan(3),r
   real, parameter :: spacing=0.05, theta_z=pi/2.
   real, parameter :: translate(3) = (/20.,0.,0./)

   nbody = 10
   narm  = 5
   nhead = 10
   nleg=5
   np = nbody + (2*narm + 1) + nhead + 2*nleg
   head_radius=spacing*narm

   allocate(xall(3,np),vall(3,np))
   vall = 0.
   xall = 0.

   index = 0

   do i=1,nbody
      index = index+1
      xall(1,index)=spacing*(i-1)
   enddo

   do i=-narm,narm
      index = index+1
      xall(:,index) = (/0.75*(spacing*nbody),spacing*i,0./)
   enddo

   do i=1,nhead
      index = index+1
      xall(:,index) = (/nbody*spacing + head_radius + head_radius*cos(i*2.*pi/10.),head_radius*sin(i*2.*pi/10.),0./)
   enddo

   do i=1,nleg
      index=index+1
      xall(1,index) = -nleg*spacing+(i-1)*spacing
      xall(2,index) = -xall(1,index)
   enddo
   do i=1,nleg
      index=index+1
      xall(1,index) = -nleg*spacing+(i-1)*spacing
      xall(2,index) = xall(1,index)
   enddo

   call get_rotation_matrix(theta_z,rotate_z,'z')

   do i=1,np
      xall(:,i)=matmul(rotate_z,xall(:,i))
      xall(:,i)=xall(:,i) + translate
      call cross_product((/0.,0.,1./),xall(:,i),rtan)
      rtan = rtan/sqrt(dot_product(rtan,rtan)) ! Unit vector tangential to motion
      r = sqrt(dot_product(xall(:,i),xall(:,i)))
      ! vall(:,i) = rtan/sqrt(r)
   enddo

   ! vall(2,:) = 0.0521157

end subroutine setup_dude

!--- Several types of single particle geodesics
subroutine initialise(x,v,type,r0)
   use metric, only: metric_type, a!,rs
   use metric_tools, only: coordinate_sys
   use force_gr, only: get_sourceterms
   use utils_gr, only: get_u0
   real, intent(in), optional  :: r0
   real, intent(out) :: x(3), v(3)
   real :: r, vy, x1
   real :: ra,va,omega,fac
   character(len=*), intent(in) :: type
   real :: rotate_y(3,3), inclination

   if (present(r0)) then
      r = r0
      write(*,'(a,f6.2)') ' Using init with r = ',r0
   endif

   select case(type)

   case('circular')
      print*,'#--- Circular velocity in x-y plane, anticlockwise ---#'
      if (.not. present(r0)) then
         print*,'Enter radius r for (r,theta,phi)=(r,pi/2,0):'
         read(*,*) r
      endif
      if (metric_type=='Schwarzschild') then
         x1 = r                                    ! x1 = r in Schwarzschild
         omega = sqrt(1./r**3)
         vy = sqrt(1./r)
      elseif(metric_type=='Kerr') then
         ! r = x1 !sqrt(x1**2-a**2)                ! x1 /= r in Kerr
         omega = 1./(r**(1.5)+a)
         x1 = sqrt(r**2 + a**2)
         vy = x1*omega
      elseif(metric_type=='Minkowski') then
         STOP 'Cannot make circular orbits in Minkowski metric.'
      endif
      select case(coordinate_sys)
      case('Cartesian')
         x = (/x1,0.,0./)
         v = (/0.,vy,0./)
      case('Spherical')
         x = (/r,0.5*pi,0./)
         v = (/0.,0.,omega/)
      end select
      print*,'period =',2*pi/omega
      print*,'Press ENTER to continue:'
      read*

   case ('radial') ! Radial infall
      print*,'#--- Radial infall ---#'
      if (.not. present(r0)) then
         print*,'Enter radius r for (r,theta,phi)=(r,pi/2,0):'
         read(*,*) r
      endif
      select case(coordinate_sys)
      case('Cartesian')
         x1 = sqrt(r**2 + a**2)
         x = (/x1,0.,0./)
      case('Spherical')
         x = (/r,0.5*pi,0./)
      end select
      v = (/0.,0.,0./)

   case('precession') ! Clement's orbit
      ra = 90.
      va = 0.0521157 ! velocity giving a pericenter rp = 10
      print*,'#--- Precessing orbit: r =',ra,' and vy = ',va
      if (present(r0)) then
         print*, "WARNING: You're trying to set r = ",r0," but r is already set for this type of orbit. Continue?"
         read*,
      endif
      select case(coordinate_sys)
      case('Cartesian')
         x = (/ra,0.,0./)
         v = (/0.,va,0./)
      case('Spherical')
         x = (/ra,0.5*pi,0./)
         v = (/0.,0.,va/ra /)
         if (.not. metric_type=='Schwarzschild') STOP 'Only have precession setup for spherical in Schwarzschild'
      end select

   case('precession inclined')
      ra = 90.
      va = 0.0521157 ! velocity giving a pericenter rp = 10
      inclination = -pi/6
      print*,'#--- Inclined precessing orbit: r =',ra,' and vy = ',va,' inclined to plane by ',inclination*180./pi, 'degrees.'
      if (present(r0)) then
         print*, "WARNING: You're trying to set r = ",r0," but r is already set for this type of orbit. Continue?"
         read*,
      endif
      if (.not. metric_type=='Kerr') then
         print*,"Warning, you are not using the Kerr metric for 'precession inclined'...result will be same as Schwarzschild"
         read*
      endif
      if (.not. coordinate_sys=='Cartesian') STOP "Haven't tested 'precession inclined' for Spherical coordinates"
      x = (/ra,0.,0./)
      v = (/0.,va,0./)
      call get_rotation_matrix(inclination,rotate_y,'y')
      x = matmul(rotate_y,x)

   case('epicycle')
      print*,'#--- Radial epicyclic motion in x-y plane, anticlockwise ---#'
      if (.not. present(r0)) then
         print*,'Enter radius r for (r,theta,phi)=(r,pi/2,0):'
         read(*,*) r
      endif
      if (metric_type=='Schwarzschild') then
         x1 = r                                    ! x1 = r in Schwarzschild
         omega = sqrt(1./r**3)
         ! vy    = sqrt(1./r)
         vy    = x1*omega
      elseif(metric_type=='Kerr') then
         ! r = x1 !sqrt(x1**2-a**2)                ! x1 /= r in Kerr
         omega = 1./(r**(1.5)+a)
         x1 = sqrt(r**2 + a**2)
         vy = x1*omega
      elseif(metric_type=='Minkowski') then
         STOP 'Cannot make circular orbits in Minkowski metric.'
      endif
      fac = 1.00001
      select case(coordinate_sys)
      case('Cartesian')
         x = (/x1,0.,0./)
         v = (/0.,fac*vy,0./)
      case('Spherical')
         x = (/r,0.5*pi,0./)
         v = (/0.,0.,fac*omega/)
      end select

   case('vertical-oscillation')
      print*,'#--- Small vertical-oscillation from circular orbit in x-y plane, anticlockwise ---#'
      if (.not. present(r0)) then
         print*,'Enter radius r for (r,theta,phi)=(r,pi/2,0):'
         read(*,*) r
      endif
      if (metric_type=='Schwarzschild') then
         x1 = r                                    ! x1 = r in Schwarzschild
         omega = sqrt(1./r**3)
         ! vy    = sqrt(1./r)
         vy    = x1*omega
      elseif(metric_type=='Kerr') then
         ! r = x1 !sqrt(x1**2-a**2)                ! x1 /= r in Kerr
         omega = 1./(r**(1.5)+a)
         x1 = sqrt(r**2 + a**2)
         vy = x1*omega
      elseif(metric_type=='Minkowski') then
         STOP 'Cannot make circular orbits in Minkowski metric.'
      endif
      fac = 1.00001
      select case(coordinate_sys)
      case('Cartesian')
         x = (/x1,0.,0.+(fac-1.)/)
         v = (/0.,vy,0./)
      case('Spherical')
         x = (/r,0.5*pi*fac,0./)
         v = (/0.,0.,omega/)
      end select

   end select
end subroutine initialise

!--- Subroutine to compute the 3-vector cross product
subroutine cross_product(a,b,c)
   real, dimension(3), intent(in) :: a,b
   real, dimension(3), intent(out) :: c

   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)

end subroutine cross_product

!--- Subroutine to return a rotation matrix, given some angle, about
!    one of the cartesian axes.
subroutine get_rotation_matrix(angle,rotation_matrix,axis)
   real, intent(in)  :: angle
   character(len=*), intent(in) :: axis
   real, intent(out) ::  rotation_matrix(3,3)
   real, dimension(3,3) :: rotate_x,rotate_y,rotate_z

   select case(axis)
   case('x')
      rotate_x(1,:) = (/1. , 0.        , 0.        /)
      rotate_x(2,:) = (/0. , cos(angle),-sin(angle)/)
      rotate_x(3,:) = (/0. , sin(angle), cos(angle)/)
      rotation_matrix = rotate_x

   case('y')
      rotate_y(1,:) = (/ cos(angle), 0. , sin(angle)/)
      rotate_y(2,:) = (/ 0.        , 1. , 0.        /)
      rotate_y(3,:) = (/-sin(angle), 0. , cos(angle)/)
      rotation_matrix = rotate_y

   case('z')
      rotate_z(1,:) = (/cos(angle),-sin(angle),0./)
      rotate_z(2,:) = (/sin(angle), cos(angle),0./)
      rotate_z(3,:) = (/0.        ,0.         ,1./)
      rotation_matrix = rotate_z
   end select

end subroutine get_rotation_matrix

end module init
