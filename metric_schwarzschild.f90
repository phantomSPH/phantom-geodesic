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
 endif

 sqrtg=1.

end subroutine get_metric

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
 real,    dimension(1+ndim,1+ndim) :: gcov, gcon
 real    :: sqrtg

 call get_metric(ndim,x,gcov,gcon,sqrtg)


end subroutine get_sourceterms

end module metric
