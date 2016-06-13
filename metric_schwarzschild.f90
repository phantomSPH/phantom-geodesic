!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2015 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: metric
!
!  DESCRIPTION:
!   Implements the Schwarzschild metric
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id: 74a92ea0d209b843dc122322eff4bee581a32acb $
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module metric
 implicit none
 character(len=*), parameter :: metric_type = 'Schwarzschild'

 public :: get_metric, get_sourceterms

 real, parameter, public :: mass1 = 1.  ! mass of central object

 private

contains

 !----------------------------------------------------------------
 !+
 !  Compute the metric tensor in both covariant (gcov) and
 !  contravariant (gcon) form
 !+
 !----------------------------------------------------------------
 pure subroutine get_metric(x,gcov,gcon,sqrtg)
  real,    intent(in)  :: x(3)
  real,    intent(out) :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
  real :: rs,r2,r,r3,rs_on_r3,coeff

  gcov = 0.
  gcon = 0.
  rs = 2.*mass1
  r2 = dot_product(x,x)
  r  = sqrt(r2)
  r3 = r**3
  rs_on_r3 = rs/r3

  if (metric_type .eq. 'Schwarzschild') then
   coeff = 1./(1.-rs/r)
   gcov(0,0) = -(1.-rs/r)
   gcov(1,1) = coeff*(1.-rs_on_r3*(x(2)**2+x(3)**2))
   gcov(2,2) = coeff*(1.-rs_on_r3*(x(1)**2+x(3)**2))
   gcov(3,3) = coeff*(1.-rs_on_r3*(x(1)**2+x(2)**2))
   gcov(1,2) = coeff*x(1)*x(2)*rs_on_r3
   gcov(2,1) = gcov(1,2)
   gcov(1,3) = coeff*x(1)*x(3)*rs_on_r3
   gcov(3,1) = gcov(1,3)
   gcov(2,3) = coeff*x(2)*x(3)*rs_on_r3
   gcov(3,2) = gcov(2,3)

   sqrtg=1.

   gcon(0,0) = -1./(1.-rs/r)
   gcon(1,1) = 1.-rs_on_r3*x(1)**2
   gcon(2,2) = 1.-rs_on_r3*x(2)**2
   gcon(3,3) = 1.-rs_on_r3*x(3)**2
   gcon(1,2) = -rs_on_r3*x(1)*x(2)
   gcon(2,1) = gcon(1,2)
   gcon(1,3) = -rs_on_r3*x(1)*x(3)
   gcon(3,1) = gcon(1,3)
   gcon(2,3) = -rs_on_r3*x(2)*x(3)
   gcon(3,2) = gcon(2,3)
  endif



 end subroutine get_metric


 !----------------------------------------------------------------
 !+
 !  Compute the derivatives of the metric tensor in both covariant
 !  (dgcovdx, dgcovdy, dgcovdz) and contravariant (dgcondx, dgcondy,
 !  dgcondz) form
 !+
 !----------------------------------------------------------------
 pure subroutine get_metric_derivs(position,dgcovdx, dgcovdy, dgcovdz, dgcondx, dgcondy, dgcondz)
  real,    intent(in)  :: position(3)
  real,    intent(out) :: dgcovdx(0:3,0:3), dgcovdy(0:3,0:3), dgcovdz(0:3,0:3)
  real,    intent(out) :: dgcondx(0:3,0:3), dgcondy(0:3,0:3), dgcondz(0:3,0:3)
  real :: x,y,z, rs,r2,r,r3,rs_on_r3

  dgcovdx = 0.
  dgcovdy = 0.
  dgcovdz = 0.
  dgcondx = 0.
  dgcondy = 0.
  dgcondz = 0.

  x = position(1)
  y = position(2)
  z = position(3)
  rs = 2.*mass1
  r2 = dot_product(position,position)
  r  = sqrt(r2)
  r3 = r**3
  rs_on_r3 = rs/r3
  if (metric_type .eq. 'Schwarzschild') then
   ! dx
   dgcovdx(0,0) = -((rs*x)/r**3)
   dgcovdx(0,1) = 0
   dgcovdx(0,2) = 0
   dgcovdx(0,3) = 0
   dgcovdx(1,0) = 0
   dgcovdx(1,1) = (x*(2*r**4 - 2*r**3*rs - 2*rs**2*(y**2 + z**2) - 2*r**2*(x**2 &
   + y**2 + z**2) + r*rs*(x**2 + 4*(y**2 + z**2))))/(r**4*(r - rs)**2)
   dgcovdx(1,2) = (rs*(r**3 - r**2*rs - 3*r*x**2 + 2*rs*x**2)*y)/(r**4*(r - rs)**2)
   dgcovdx(1,3) = (rs*(r**3 - r**2*rs - 3*r*x**2 + 2*rs*x**2)*z)/(r**4*(r - rs)**2)
   dgcovdx(2,0) = 0
   dgcovdx(2,1) = (rs*(r**3 - r**2*rs - 3*r*x**2 + 2*rs*x**2)*y)/(r**4*(r - rs)**2)
   dgcovdx(2,2) = (x*(2*r**4 - 4*r**3*rs + 2*r**2*(rs**2 - x**2 - y**2 - z**2) &
   - 2*rs**2*(x**2 + z**2) + r*rs*(4*x**2 + y**2 + 4*z**2)))/(r**4*(r - rs)**2)
   dgcovdx(2,3) = (rs*(-3*r + 2*rs)*x*y*z)/(r**4*(r - rs)**2)
   dgcovdx(3,0) = 0
   dgcovdx(3,1) = (rs*(r**3 - r**2*rs - 3*r*x**2 + 2*rs*x**2)*z)/(r**4*(r - rs)**2)
   dgcovdx(3,2) = (rs*(-3*r + 2*rs)*x*y*z)/(r**4*(r - rs)**2)
   dgcovdx(3,3) = (-2*(r - rs)**2*x*(-r**2 + x**2 + y**2) + r*(-2*r + rs)*x*z**2)&
   /(r**4*(r - rs)**2)

   dgcondx(0,0) = (rs*x)/(r*(r - rs)**2)
   dgcondx(0,1) = 0
   dgcondx(0,2) = 0
   dgcondx(0,3) = 0
   dgcondx(1,0) = 0
   dgcondx(1,1) = -((x*(-2*r**3 + 2*r**2*rs - 3*rs*x**2 + 2*r*(x**2 + y**2 + z**2)))/r**5)
   dgcondx(1,2) = -((rs*(r**2 - 3*x**2)*y)/r**5)
   dgcondx(1,3) = -((rs*(r**2 - 3*x**2)*z)/r**5)
   dgcondx(2,0) = 0
   dgcondx(2,1) = -((rs*(r**2 - 3*x**2)*y)/r**5)
   dgcondx(2,2) = (x*(2*r**3 + 3*rs*y**2 - 2*r*(x**2 + y**2 + z**2)))/r**5
   dgcondx(2,3) = (3*rs*x*y*z)/r**5
   dgcondx(3,0) = 0
   dgcondx(3,1) = -((rs*(r**2 - 3*x**2)*z)/r**5)
   dgcondx(3,2) = (3*rs*x*y*z)/r**5
   dgcondx(3,3) = (x*(2*r**3 + 3*rs*z**2 - 2*r*(x**2 + y**2 + z**2)))/r**5

   ! dy
   dgcovdy(0,0) = -((rs*y)/r**3)
   dgcovdy(0,1) = 0
   dgcovdy(0,2) = 0
   dgcovdy(0,3) = 0
   dgcovdy(1,0) = 0
   dgcovdy(1,1) = (y*(2*r**4 - 4*r**3*rs + 2*r**2*(rs**2 - x**2 - y**2 - z**2) &
   - 2*rs**2*(y**2 + z**2) + r*rs*(x**2 + 4*(y**2 + z**2))))/(r**4*(r - rs)**2)
   dgcovdy(1,2) = (rs*x*(r**3 - r**2*rs - 3*r*y**2 + 2*rs*y**2))/(r**4*(r - rs)**2)
   dgcovdy(1,3) = (rs*(-3*r + 2*rs)*x*y*z)/(r**4*(r - rs)**2)
   dgcovdy(2,0) = 0
   dgcovdy(2,1) = (rs*x*(r**3 - r**2*rs - 3*r*y**2 + 2*rs*y**2))/(r**4*(r - rs)**2)
   dgcovdy(2,2) = (y*(2*r**4 - 2*r**3*rs - 2*rs**2*(x**2 + z**2) - 2*r**2*(x**2&
   + y**2 + z**2) + r*rs*(4*x**2 + y**2 + 4*z**2)))/(r**4*(r - rs)**2)
   dgcovdy(2,3) = (rs*(r**3 - r**2*rs - 3*r*y**2 + 2*rs*y**2)*z)/(r**4*(r - rs)**2)
   dgcovdy(3,0) = 0
   dgcovdy(3,1) = (rs*(-3*r + 2*rs)*x*y*z)/(r**4*(r - rs)**2)
   dgcovdy(3,2) = (rs*(r**3 - r**2*rs - 3*r*y**2 + 2*rs*y**2)*z)/(r**4*(r - rs)**2)
   dgcovdy(3,3) = (-2*(r - rs)**2*y*(-r**2 + x**2 + y**2) + r*(-2*r + rs)*y*z**2)&
   /(r**4*(r - rs)**2)

   dgcondy(0,0) = (rs*y)/(r*(r - rs)**2)
   dgcondy(0,1) = 0
   dgcondy(0,2) = 0
   dgcondy(0,3) = 0
   dgcondy(1,0) = 0
   dgcondy(1,1) = (y*(2*r**3 + 3*rs*x**2 - 2*r*(x**2 + y**2 + z**2)))/r**5
   dgcondy(1,2) = -((rs*x*(r**2 - 3*y**2))/r**5)
   dgcondy(1,3) = (3*rs*x*y*z)/r**5
   dgcondy(2,0) = 0
   dgcondy(2,1) = -((rs*x*(r**2 - 3*y**2))/r**5)
   dgcondy(2,2) = -((y*(-2*r**3 + 2*r**2*rs - 3*rs*y**2 + 2*r*(x**2 + y**2 + z**2)))/r**5)
   dgcondy(2,3) = -((rs*(r**2 - 3*y**2)*z)/r**5)
   dgcondy(3,0) = 0
   dgcondy(3,1) = (3*rs*x*y*z)/r**5
   dgcondy(3,2) = -((rs*(r**2 - 3*y**2)*z)/r**5)
   dgcondy(3,3) = (y*(2*r**3 + 3*rs*z**2 - 2*r*(x**2 + y**2 + z**2)))/r**5

   ! dz
   dgcovdz(0,0) = -((rs*z)/r**3)
   dgcovdz(0,1) = 0
   dgcovdz(0,2) = 0
   dgcovdz(0,3) = 0
   dgcovdz(1,0) = 0
   dgcovdz(1,1) = (z*(2*r**4 - 4*r**3*rs + 2*r**2*(rs**2 - x**2 - y**2 - z**2) &
   - 2*rs**2*(y**2 + z**2) + r*rs*(x**2 + 4*(y**2 + z**2))))/(r**4*(r - rs)**2)
   dgcovdz(1,2) = (rs*(-3*r + 2*rs)*x*y*z)/(r**4*(r - rs)**2)
   dgcovdz(1,3) = (rs*x*(r**3 - r**2*rs - 3*r*z**2 + 2*rs*z**2))/(r**4*(r - rs)**2)
   dgcovdz(2,0) = 0
   dgcovdz(2,1) = (rs*(-3*r + 2*rs)*x*y*z)/(r**4*(r - rs)**2)
   dgcovdz(2,2) = (z*(2*r**4 - 4*r**3*rs + 2*r**2*(rs**2 - x**2 - y**2 - z**2) -&
   2*rs**2*(x**2 + z**2) + r*rs*(4*x**2 + y**2 + 4*z**2)))/(r**4*(r - rs)**2)
   dgcovdz(2,3) = (rs*y*(r**3 - r**2*rs - 3*r*z**2 + 2*rs*z**2))/(r**4*(r - rs)**2)
   dgcovdz(3,0) = 0
   dgcovdz(3,1) = (rs*x*(r**3 - r**2*rs - 3*r*z**2 + 2*rs*z**2))/(r**4*(r - rs)**2)
   dgcovdz(3,2) = (rs*y*(r**3 - r**2*rs - 3*r*z**2 + 2*rs*z**2))/(r**4*(r - rs)**2)
   dgcovdz(3,3) = (z*(2*(r - rs)*(r**3 - r*(x**2 + y**2) + rs*(x**2 + y**2)) + &
   r*(-2*r + rs)*z**2))/(r**4*(r - rs)**2)

   dgcondz(0,0) = (rs*z)/(r*(r - rs)**2)
   dgcondz(0,1) = 0
   dgcondz(0,2) = 0
   dgcondz(0,3) = 0
   dgcondz(1,0) = 0
   dgcondz(1,1) = (z*(2*r**3 + 3*rs*x**2 - 2*r*(x**2 + y**2 + z**2)))/r**5
   dgcondz(1,2) = (3*rs*x*y*z)/r**5
   dgcondz(1,3) = -((rs*x*(r**2 - 3*z**2))/r**5)
   dgcondz(2,0) = 0
   dgcondz(2,1) = (3*rs*x*y*z)/r**5
   dgcondz(2,2) = (z*(2*r**3 + 3*rs*y**2 - 2*r*(x**2 + y**2 + z**2)))/r**5
   dgcondz(2,3) = -((rs*y*(r**2 - 3*z**2))/r**5)
   dgcondz(3,0) = 0
   dgcondz(3,1) = -((rs*x*(r**2 - 3*z**2))/r**5)
   dgcondz(3,2) = -((rs*y*(r**2 - 3*z**2))/r**5)
   dgcondz(3,3) = (z*(2*r**3 - 2*r**2*rs + 3*rs*z**2 - 2*r*(x**2 + y**2 + z**2)))/r**5

  endif

 end subroutine get_metric_derivs

 !----------------------------------------------------------------
 !+
 !  Compute the source terms required on the right hand side of
 !  the relativistic momentum equation. These are of the form:
 !   T^\mu\nu dg_\mu\nu/dx^i
 !+
 !----------------------------------------------------------------
 pure subroutine get_sourceterms(x,v,fterm)
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

end module metric
