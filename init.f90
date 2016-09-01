module init
   implicit none
   real, parameter :: pi = acos(-1.)
contains

   subroutine setup(xall,vall,np)
      integer, intent(out) :: np
      real, allocatable, intent(inout), dimension(:,:) :: xall,vall
      call setup_dude(xall,vall,np)

   end subroutine setup

   subroutine setup_2testbodies(xall,vall,np)
      integer, intent(out) :: np
      real, allocatable, intent(inout), dimension(:,:) :: xall,vall
      np = 2
      allocate(xall(3,np),vall(3,np))

      call initialise(xall(:,1),vall(:,1),'precession')
      call initialise(xall(:,2),vall(:,2),'circular')

   end subroutine setup_2testbodies

   subroutine setup_dude(xall,vall,np)
      integer, intent(out) :: np
      real, allocatable, intent(inout), dimension(:,:) :: xall,vall
      integer :: i,index,nleg,nbody,narm,nhead
      real :: rotate_z(3,3),head_radius,rtan(3),r
      real, parameter :: spacing=0.05, theta_z=pi/2.
      real, parameter :: translate(3) = (/10.,0.,0./)

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

      rotate_z(1,:)=(/cos(theta_z),-sin(theta_z),0./)
      rotate_z(2,:)=(/sin(theta_z), cos(theta_z),0./)
      rotate_z(3,:)=(/0.        ,0.         ,1./)

      do i=1,np
         xall(:,i)=matmul(rotate_z,xall(:,i))
         xall(:,i)=xall(:,i) + translate
         call cross_product((/0.,0.,1./),xall(:,i),rtan)
         rtan = rtan/sqrt(dot_product(rtan,rtan))
         r = sqrt(dot_product(xall(:,i),xall(:,i)))
         vall(:,i) = rtan/sqrt(r)
      enddo

      ! vall(2,:) = 0.2

   end subroutine setup_dude

   subroutine initialise(x,v,type)
      real, intent(out) :: x(3), v(3)
      real :: r, vtan
      real :: ra,va
      ! real, parameter :: pi=3.14159265358979

      !character(len=*), parameter :: type = 'precession'
      character(len=*), intent(in) :: type

      if (type=='circular') then
         print*,'(Circular velocity) Enter radius r:'
         read(*,*) r
         vtan = sqrt(1./r)
         x = (/r,0.,0./)
         v = (/0.,vtan,0./)
         print*,'period =',2*pi*r/vtan
         print*,'Press ENTER to continue:'
         read(*,*)
      endif

      if (type=='radial') then
         ! Radial infall
         x = (/-20.,0.,0./)
         v = (/0.,0.,0./)
      endif

      if (type=='precession') then
         ! Clement's orbit
         ra = 90.
         va = 0.0521157 ! velocity giving a pericenter rp = 10
         x = (/ra,0.,0./)
         v = (/0.,va,0./)
      endif

      ! x = (/84.6,0.,0./)
      ! v = (/0.,0.04777816206847369,0./)

      ! x = (/500.,0.,0./)
      ! v = (/0.,0.044721359549996,0./)

      ! x = (/6.001,0.,0./)
      ! v = (/0.,0.136083276348796,0./)

      ! x = (/3.,0.,0./)
      ! v = (/0.,0.19245008972987526,0./)

      ! !!Circular orbit....?
      ! x = (/10.,0.,0./)
      ! v = (/0.,0.316228,0./)

   end subroutine initialise

   subroutine setup_sphere(xall,vall,np)
      integer, intent(out) :: np
      real, allocatable, intent(inout), dimension(:,:) :: xall,vall
      real :: dlayer, dr
      integer :: nr, ntheta, nphi, i,j,k, nlayers, nrings,index,n,nringsmax
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
      do i=1,np
         xall(:,i)=xall(:,i) + translate
         vall(:,i)=(/0.,0.0521157,0./)
      enddo
   end subroutine setup_sphere

   subroutine cross_product(a,b,c)
      real, dimension(3), intent(in) :: a,b
      real, dimension(3), intent(out) :: c

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)

   end subroutine cross_product

end module init
