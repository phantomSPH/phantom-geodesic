module init
 implicit none
contains
 subroutine initialise(x,v)
  real, intent(out) :: x(3), v(3)
  real :: r, vtan
  real :: ra,va

  ! print*,'read r:'
  ! read(*,*) r
  ! vtan = sqrt(1./r)
  ! x = (/r,0.,0./)
  ! v = (/0.,vtan,0./)

  ! x = (/84.6,0.,0./)
  ! v = (/0.,0.04777816206847369,0./)



  ra = 90.
  va = 0.0521157 ! velocity giving a pericenter rp = 10

  x = (/ra,0.,0./)
  v = (/0.,va,0./)

  ! x = (/84.6,0.,0./)
  ! v = (/0.,0.04777816206847369,0./)

  ! x = (/500.,0.,0./)
  ! v = (/0.,0.044721359549996,0./)

  ! x = (/6.001,0.,0./)
  ! v = (/0.,0.136083276348796,0./)

  ! x = (/3.,0.,0./)
  ! v = (/0.,0.19245008972987526,0./)

 end subroutine initialise

end module init
