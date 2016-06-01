module init
 implicit none
contains
 subroutine initialise(x,v)
  real(8), intent(out) :: x(3), v(3)

  x = (/10.,0.,0./)
  v = (/0.,0.,0./)
 end subroutine initialise

end module init
