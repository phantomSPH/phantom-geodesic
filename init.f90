module init
   implicit none
contains

   subroutine setup(xall,vall,np)
      integer, intent(out) :: np
      real, allocatable, intent(inout), dimension(:,:) :: xall,vall
      np = 2
      allocate(xall(3,np),vall(3,np))

      call initialise(xall(:,1),vall(:,1),'precession')
      call initialise(xall(:,2),vall(:,2),'circular')

   end subroutine setup

   subroutine setup_dude(xall,vall,np)
      integer, intent(out) :: np
      real, allocatable, intent(inout), dimension(:,:) :: xall,vall
      integer :: i,index,nleg,nbody,narm,nhead
      real, parameter :: spacing=1.,head_radius=10., pi = acos(-1.), start=90.

      nbody = 10
      narm  = 5
      nhead = 10
      nleg=5
      np = nbody + (2*narm + 1) + nhead + 2*nleg

      allocate(xall(3,np),vall(3,np))
      vall = 0.
      xall = 0.

      index = 0

      do i=1,nbody
         index = index+1
         xall(1,index)=start+spacing*(i-1)
      enddo

      do i=-narm,narm
         index = index+1
         xall(:,index) = (/0.75*(spacing*nbody) + start,spacing*i,0./)
      enddo

      do i=1,nhead
         index = index+1
         xall(:,index) = (/start+nbody*spacing + head_radius + head_radius*cos(i*2.*pi/10.),head_radius*sin(i*2.*pi/10.),0./)
      enddo

      do i=1,nleg
         index=index+1
         xall(1,index) = start-nleg*spacing+(i-1)*spacing
         xall(2,index) = -xall(1,index) + start
      enddo
      do i=1,nleg
         index=index+1
         xall(1,index) = start-nleg*spacing+(i-1)*spacing
         xall(2,index) = xall(1,index) - start
      enddo

      vall(2,:) = 0.0521157

   end subroutine setup_dude

   subroutine initialise(x,v,type)
      real, intent(out) :: x(3), v(3)
      real :: r, vtan
      real :: ra,va
      real, parameter :: pi=3.14159265358979

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

end module init
