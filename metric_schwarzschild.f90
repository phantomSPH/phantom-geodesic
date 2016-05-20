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
pure subroutine get_metric(ndim,x,gcov,gcon,sqrtg)
 integer, intent(in)  :: ndim
 real,    intent(in)  :: x(ndim)
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
    gcov(1,1) = coeff*(1.-rs/r3*(x(2)**2+x(3)**2))
    gcov(2,2) = coeff*(1.-rs/r3*(x(1)**2+x(3)**2))
    gcov(3,3) = coeff*(1.-rs/r3*(x(1)**2+x(2)**2))
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
    gcon(2,1) = gcov(1,2)
    gcon(1,3) = -rs_on_r3*x(1)*x(3)
    gcon(3,1) = gcov(1,3)
    gcon(2,3) = -rs_on_r3*x(2)*x(3)
    gcon(3,2) = gcov(2,3)
 endif



end subroutine get_metric


!----------------------------------------------------------------
!+
!  Compute the derivatives of the metric tensor in both covariant
!  (dgcovdx, dgcovdy, dgcovdz) and contravariant (dgcondx, dgcondy,
!  dgcondz) form
!+
!----------------------------------------------------------------
pure subroutine get_metric_derivs(ndim,x,dgcovdx, dgcovdy, dgcovdz, dgcondx, dgcondy, dgcondz)
 integer, intent(in)  :: ndim
 real,    intent(in)  :: x(ndim)
 real,    intent(out) :: dgcovdx(0:3,0:3), dgcovdy(0:3,0:3), dgcovdz(0:3,0:3)
 real,    intent(out) :: dgcondx(0:3,0:3), dgcondy(0:3,0:3), dgcondz(0:3,0:3)


end subroutine get_metric_derivs

!----------------------------------------------------------------
!+
!  Compute the source terms required on the right hand side of
!  the relativistic momentum equation. These are of the form:
!   T^\mu\nu dg_\mu\nu/dx^i
!+
!----------------------------------------------------------------
pure subroutine get_sourceterms(ndim,x,v,fterm)
 integer, intent(in)  :: ndim
 real,    intent(in)  :: x(ndim),v(ndim)
 real,    intent(out) :: fterm(ndim)
 real,    :: gcov(0:3,0:3), gcon(0:3,0:3)
 real,    :: sqrtg
 real,    :: dgcovdx(0:3,0:3), dgcovdy(0:3,0:3), dgcovdz(0:3,0:3)
 real,    :: dgcondx(0:3,0:3), dgcondy(0:3,0:3), dgcondz(0:3,0:3)
 real,    :: v4(0:3), tmunu(0:3,0:3)
 real,    :: P, rho, u, hcur, uzero2, uzero

 call get_metric(ndim,x,gcov,gcon,sqrtg)
 call get_metric_derivs(ndim,x,dgcovdx, dgcovdy, dgcovdz, dgcondx, dgcondy, dgcondz)

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
 do i=0,3
    do j=0,3
       tmunu(i,j) = rho*hcur*uzero2*v4(i)*v4(j) + P*gcont(i,j)
    enddo
 enddo

 ! source term
 fterm(1) = 0.
 fterm(2) = 0.
 fterm(3) = 0.
 do i=0,3
    do j=0,3
       fterm(1) = fterm(1) + tmunu(i,j)*dgcovdx(i,j)
       fterm(2) = fterm(2) + tmunu(i,j)*dgcovdy(i,j)
       fterm(3) = fterm(3) + tmunu(i,j)*dgcovdz(i,j)
    enddo
 enddo

 do i=1,ndim
    fterm(i) = 0.5*fterm(i) / (rho*uzero)
 enddo

end subroutine get_sourceterms

end module metric
