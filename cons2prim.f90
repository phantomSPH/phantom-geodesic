!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2015 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: cons2prim
!
!  DESCRIPTION:
!   Solve for the primitive variables from the conserved variables
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
module cons2prim
 implicit none

 public :: conservative2primitive, primitive2conservative
 
 integer, parameter :: ndim = 3

 private

contains

!----------------------------------------------------------------
!+
!  Construct conserved variables from the primitive variables
!  primitive variables are (v^i,d,u,P); v i
!  conserved variables are (rho,pmom_i,en)
!+
!----------------------------------------------------------------
subroutine primitive2conservative(x,den,v,u,P,rho,pmom,en)
! use metric,   only:get_metric
! use utils_gr, only:dot_product_gr
 real, intent(in)  :: x(ndim)
 real, intent(in)  :: den,v(ndim),u,P
 real, intent(out) :: rho,pmom(ndim),en

! call get_metric(ndim,x,gcov,gcon,sqrtg)

! Lorentz = 1./sqrt(1. + dot_product_gr(v,v,gcov))
! w = 1. + u + P/den

! rho = d*gamma*sqrtg
! pmom = v*gamma

 return
end subroutine primitive2conservative

subroutine conservative2primitive(x,v,den,u,P,rho,pmom,en)
 real, intent(in)  :: x(ndim)
 real, intent(out) :: v(ndim),den,u,P
 real, intent(in)  :: rho,pmom(ndim),en
! use metric, only:get_metric

! call get_metric(ndim,x,gcov,gcon,sqrtg)

! Lorentz = 1.
! w = 1. + u + P/den

! rho = d*gamma*sqrtg
! pmom = v*gamma

 return
end subroutine conservative2primitive

end module cons2prim
