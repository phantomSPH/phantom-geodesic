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
subroutine get_sourceterms(x,v,fterm)
 use metric, only: get_metric, get_metric_derivs
 use utils_gr, only: dot_product_gr
 real,    intent(in)  :: x(3),v(3)
 real,    intent(out) :: fterm(3)
 real    :: gcov(0:3,0:3), gcon(0:3,0:3)
 real    :: sqrtg
 real    :: dgcovdx(0:3,0:3), dgcovdy(0:3,0:3), dgcovdz(0:3,0:3)
 real    :: dgcondx(0:3,0:3), dgcondy(0:3,0:3), dgcondz(0:3,0:3)
 real    :: v4(0:3), tmunu(0:3,0:3)
 real    :: P, rho, u, hcur, uzero2, uzero, r, r2
 integer :: i,j
 ! Note to self: try with potential from Tejeda, Rosswog  2013
 call get_metric(x,gcov,gcon,sqrtg)
 call get_metric_derivs(x,dgcovdx, dgcovdy, dgcovdz, dgcondx, dgcondy, dgcondz)

 ! for a test particle
 P = 0.
 u = 0.
 rho = 1. ! this value does not matter (will cancel in the momentum equation)
 hcur = 1.

 !hcur = 1 + u + P/rho

 ! lower-case 4-velocity
 do i=1,3
  v4(i) = v(i)
 enddo
 v4(0) = 1.

 ! first component of the upper-case 4-velocity
 uzero2 = -1./dot_product_gr(v4,v4,gcov)
 uzero = sqrt(uzero2)

 ! energy-momentum tensor
 do j=0,3
  do i=0,3
   tmunu(i,j) = rho*hcur*uzero2*v4(i)*v4(j) + P*gcon(i,j)
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

 r2 = dot_product(x,x)
 r  = sqrt(r2)
 do i=1,3
  fterm(i) = 0.5*fterm(i) / (rho*uzero)
  !fterm(i) = -x(i)/(r2*r) !Keplerian?
 enddo

end subroutine get_sourceterms

end module force_gr