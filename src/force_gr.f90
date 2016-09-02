module force_gr
   implicit none
contains

!----------------------------------------------------------------
!+
!  Compute the source terms required on the right hand side of
!  the relativistic momentum equation. These are of the form:
!   T^\mu\nu dg_\mu\nu/dx^i
!+
!----------------------------------------------------------------
subroutine get_forcegr(x,v,dens,u,p,fterm)
 use metric, only: get_metric, get_metric_derivs
 use utils_gr, only: get_u0
 use eos, only: get_enthalpy
 real,    intent(in)  :: x(3),v(3),dens,u,p
 real,    intent(out) :: fterm(3)
 real    :: gcov(0:3,0:3), gcon(0:3,0:3)
 real    :: sqrtg
 real    :: dgcovdx(0:3,0:3), dgcovdy(0:3,0:3), dgcovdz(0:3,0:3)
 ! real    :: dgcondx(0:3,0:3), dgcondy(0:3,0:3), dgcondz(0:3,0:3)
 real    :: v4(0:3), tmunu(0:3,0:3)
 real    :: enth, uzero2, uzero
 integer :: i,j
 ! Note to self: try with potential from Tejeda, Rosswog  2013
 call get_metric(x,gcov,gcon,sqrtg)
 call get_metric_derivs(x,dgcovdx, dgcovdy, dgcovdz)

 call get_enthalpy(enth,dens,p)

 ! lower-case 4-velocity
 v4(0) = 1.
 v4(1:3) = v(:)

 ! first component of the upper-case 4-velocity
 call get_u0(x,v,uzero)
 uzero2 = uzero**2

 ! energy-momentum tensor
 do j=0,3
  do i=0,3
   tmunu(i,j) = dens*enth*uzero2*v4(i)*v4(j) + P*gcon(i,j)
  enddo
 enddo

 ! source term
 fterm(1) = 0.
 fterm(2) = 0.
 fterm(3) = 0.
 do j=0,3
  do i=0,3
   fterm(1) = fterm(1) + tmunu(i,j)*dgcovdx(i,j)
   fterm(2) = fterm(2) + tmunu(i,j)*dgcovdy(i,j)
   fterm(3) = fterm(3) + tmunu(i,j)*dgcovdz(i,j)
  enddo
 enddo

 fterm = 0.5*fterm/(dens*uzero)

end subroutine get_forcegr

! for a test particle
subroutine get_sourceterms(x,v,fterm)
 real, intent(in)  :: x(3),v(3)
 real, intent(out) :: fterm(3)
 real :: dens,u,p

 P = 0.
 u = 0.
 dens = 1. ! this value does not matter (will cancel in the momentum equation)
 call get_forcegr(x,v,dens,u,p,fterm)
end subroutine get_sourceterms

end module force_gr
