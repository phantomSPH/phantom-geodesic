module forces
 implicit none
contains

!----------------------------------------------------------------
!+
!
!  calculate gravitational force between particles
!
!+
!----------------------------------------------------------------
subroutine get_newtonian_force_new(np,xall,fterm_new,mall,i)
  real, dimension(:,:), intent(in) :: xall
  real, dimension(:), intent(in) :: mall
  integer, intent(in) :: np,i
  real, dimension(3), intent(inout) :: fterm_new
  real, dimension(3) :: x,x_other,x_rel
  real :: rr,ddr
  integer :: j

    ! fterm_new(:,i) = 0.
    ! we determine the force on each particle
    x = xall(:,i)
    do j = 1,np
      if (i .ne. j) then
        ! determine relative position
        x_other = xall(:,j)
        x_rel = x - x_other
        rr = sqrt(dot_product(x_rel,x_rel))
        ddr = 1.0 / rr**3
        ! determine the force components in vector form
        fterm_new(:) = fterm_new(:) - ddr*x_rel*mall(j)
      endif
    enddo
    ! print*, "fterm in the forces file: ", fterm_new

end subroutine get_newtonian_force_new

end module forces
