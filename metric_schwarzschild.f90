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
 real,    intent(out) :: gcov(ndim,ndim), gcon(ndim,ndim), sqrtg

end subroutine get_metric

!----------------------------------------------------------------
!+
!  Compute the source terms required on the right hand side of
!  the relativistic momentum equation. These are of the form:
!   T^\mu\nu dg_\mu\nu/dx^i
!+
!----------------------------------------------------------------
pure subroutine get_sourceterms(ndim,x,v,gcov,gcon,fterm)
 integer, intent(in)  :: ndim
 real,    intent(in)  :: x(ndim),v(ndim)
 real,    intent(in)  :: gcov(ndim,ndim), gcon(ndim,ndim)
 real,    intent(out) :: fterm(ndim)
 
end subroutine get_sourceterms

end module metric
