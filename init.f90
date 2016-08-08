module init
 implicit none
contains
 subroutine initialise(x,v)
  real, intent(out) :: x(3), v(3)
  real :: r, vtan
  real :: ra,va
  real, parameter :: pi=3.14159265358979
  character(len=*), parameter :: type = 'circular'

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
