!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2015 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: utils_gr
!
!  DESCRIPTION:
!   Contains utility routines for general relativity
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
module utils_gr
 implicit none

 public :: dot_product_gr, get_metric3plus1, get_u0

 private

contains

!----------------------------------------------------------------
!+
!  Function to perform a dot product in general relativity
!  i.e. g_\mu\nu v^\mu \v^nu
!+
!----------------------------------------------------------------
pure real function dot_product_gr(vec1,vec2,gcov)
 real, intent(in) :: vec1(:)
 real, intent(in) :: vec2(size(vec1))
 real, intent(in) :: gcov(size(vec1),size(vec2))
 integer :: i,j

 dot_product_gr = 0.
 do i=1,size(vec1)
    do j=1,size(vec2)
       dot_product_gr = dot_product_gr + gcov(j,i)*vec1(i)*vec2(j)
    enddo
 enddo

 return
end function dot_product_gr

subroutine get_metric3plus1(x,alpha,beta,gij)
   use metric, only: get_metric
   real, intent(in) :: x(1:3)
   real, intent(out) :: alpha,beta(1:3),gij(1:3,1:3)
   real :: gcov(0:3,0:3),gcon(0:3,0:3),sqrtg,beta2

   call get_metric(x,gcov,gcon,sqrtg)
   beta  = gcov(0,1:3)
   beta2 = dot_product_gr(beta,beta,gcon(1:3,1:3))
   alpha = sqrt(beta2 - gcov(0,0))
   gij   = gcov(1:3,1:3)

end subroutine get_metric3plus1

subroutine get_u0(x,v,U0)
   use metric, only: get_metric
   real, intent(in) :: x(1:3),v(1:3)
   real, intent(out) :: U0
   real :: v4(0:3), gcov(0:3,0:3), gcon(0:3,0:3), sqrtg

   call get_metric(x,gcov,gcon,sqrtg)
   v4(0) = 1.
   v4(1:3) = v(1:3)
   U0 = 1./sqrt(-dot_product_gr(v4,v4,gcov))
end subroutine get_u0
end module utils_gr
