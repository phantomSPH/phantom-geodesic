module init
 implicit none
contains
 subroutine initialise(x,v)
  real(8), intent(out) :: x(3), v(3)

  !x = (/1.d0,0.d0,0.d0/)
  !v = (/0.d0,1.d0,0.d0/)
 end subroutine initialise

end module init
