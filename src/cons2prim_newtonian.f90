module cons2prim
 implicit none

 public :: primitive2conservative, conservative2primitive,get_p_from_v, get_v_from_p
 public :: error_to_string
 integer, parameter :: ierr_notconverged = 1

 private

contains

elemental function error_to_string(ierr) result(string)
 integer, intent(in) :: ierr
 character(len=20) :: string

 select case(ierr)
 case(ierr_notconverged)
    string = 'not converged'
 case default
    string = 'unknown error'
 end select
end function error_to_string

!----------------------------------------------------------------
!+
!  Construct conserved variables from the primitive variables
!  primitive variables are (v^i,d,u,P); v i
!  conserved variables are (rho,pmom_i,en)
!+
!----------------------------------------------------------------
subroutine primitive2conservative(x,v,dens,u,P,rho,pmom,en,en_type)
 real, intent(in)  :: x(1:3)
 real, intent(in) :: dens,v(1:3),u,P
 real, intent(out)  :: rho,pmom(1:3),en
 character(len=*), intent(in) :: en_type

end subroutine primitive2conservative

subroutine conservative2primitive(x,v,dens,u,P,rho,pmom,en,ierr,en_type)
 real, intent(in)  :: x(1:3)
 real, intent(inout) :: v(1:3),dens,u,P
 real, intent(in)  :: rho,pmom(1:3),en
 integer, intent(out) :: ierr
 character(len=*), intent(in) :: en_type


end subroutine conservative2primitive

subroutine get_v_from_p(pmom,v,x)
 real, intent(in) :: pmom(1:3), x(1:3)
 real, intent(out) :: v(1:3)

 v = pmom

end subroutine get_v_from_p

subroutine get_p_from_v(pmom,v,x)
 real, intent(in) :: v(1:3), x(1:3)
 real, intent(out) :: pmom(1:3)

 pmom = v

end subroutine get_p_from_v

end module cons2prim
