module init
 implicit none
contains
 subroutine initialise(ndim,x,v)
  integer, intent(in) :: ndim
  real(8), intent(out) :: x(ndim), v(ndim)

  !x = (/1.d0,0.d0,0.d0/)
  !v = (/0.d0,1.d0,0.d0/)
 end subroutine initialise

end module init
