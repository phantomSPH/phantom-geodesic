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

 public :: conservative2primitive, primitive2conservative, get_p_from_v, get_v_from_p

 integer, parameter :: ndim = 4

 private

contains

!----------------------------------------------------------------
!+
!  Construct conserved variables from the primitive variables
!  primitive variables are (v^i,d,u,P); v i
!  conserved variables are (rho,pmom_i,en)
!+
!----------------------------------------------------------------
subroutine primitive2conservative(x,v,den,u,P,rho,pmom,en)
 use utils_gr, only: dot_product_gr
 use metric, only: get_metric
 real, intent(in)  :: x(1:3)
 real, intent(in) :: den,v(1:3),u,P
 real, intent(out)  :: rho,pmom(1:3),en
 real, dimension(0:3,0:3) :: gcov, gcon
 real :: sqrtg, H, gimuvmu, gimuvmuvi, U0, v4(0:3)
 integer :: i, mu
 
 v4 = (/1.,v(1:3)/)
 H = 1.+ u + P/den
 call get_metric(x,gcov,gcon,sqrtg)
 U0  = 1./sqrt(-dot_product_gr(v,v,gcov))

 rho = sqrtg*den*U0
 do i=1,3
    gimuvmu = dot_product(gcov(i,:),v4(:))
    pmom(i) = U0*H*gimuvmu
 enddo
 gimuvmuvi = 0.
 do i=1,3
    do mu=1,4
       gimuvmuvi = gimuvmuvi + gcov(i,mu)*v(mu)*v(i)
    enddo
 enddo
 en = U0*(H*gimuvmuvi - (1.+u)*dot_product_gr(v,v,gcov))

 return
end subroutine primitive2conservative

subroutine conservative2primitive(x,v,den,u,P,rho,pmom,en)
 real, intent(in)  :: x(1:3)
 real, intent(out) :: v(1:3),den,u,P
 real, intent(in)  :: rho,pmom(1:3),en
! use metric, only:get_metric

! call get_metric(ndim,x,gcov,gcon,sqrtg)

! Lorentz = 1.
! w = 1. + u + P/den

! rho = d*gamma*sqrtg
! pmom = v*gamma

 return
end subroutine conservative2primitive

subroutine get_v_from_p(pmom,v,x)
 use metric, only: get_metric
 real, intent(in) :: pmom(1:3), x(1:3)
 real, intent(out) :: v(1:3)
 real :: den, u, en, P, rho

 en = 0. ! ???
 P  = 0.
 rho = 0. ! ???
 call conservative2primitive(x,v,den,u,P,rho,pmom,en)

 return
end subroutine get_v_from_p

subroutine get_p_from_v(pmom,v,x)
 use metric, only: get_metric
 real, intent(in) :: v(1:3), x(1:3)
 real, intent(out) :: pmom(1:3)
 real :: rho, en, den, u, P

 den = 1. ! ???
 u = 0.
 P = 0.
 call primitive2conservative(x,v,den,u,P,rho,pmom,en)

 return
end subroutine get_p_from_v

end module cons2prim
