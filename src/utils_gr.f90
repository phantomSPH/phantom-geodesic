module utils_gr
 implicit none

 public :: dot_product_gr, get_metric3plus1, get_u0, get_rderivs

 private
 interface get_metric3plus1
  module procedure get_metric3plus1_only, get_metric3plus1_both
 end interface get_metric3plus1

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
 real :: vec1i
 integer :: i,j

 dot_product_gr = 0.
 do i=1,size(vec1)
    vec1i = vec1(i)
    do j=1,size(vec2)
       dot_product_gr = dot_product_gr + gcov(j,i)*vec1i*vec2(j)
    enddo
 enddo

 return
end function dot_product_gr

subroutine get_metric3plus1_only(x,alpha,beta,gammaijdown,gammaijUP)
 real, intent(in)  :: x(1:3)
 real, intent(out) :: alpha, beta(1:3), gammaijdown(1:3,1:3),gammaijUP(1:3,1:3)
 real              :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg

 call metric3p1(x,alpha,beta,gammaijdown,gammaijUP,gcov,gcon,sqrtg)

end subroutine get_metric3plus1_only

subroutine get_metric3plus1_both(x,alpha,beta,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
 real, intent(in)  :: x(1:3)
 real, intent(out) :: alpha,beta(1:3), gammaijdown(1:3,1:3),gammaijUP(1:3,1:3)
 real, intent(out) :: gcov(0:3,0:3),gcon(0:3,0:3),sqrtg

 call metric3p1(x,alpha,beta,gammaijdown,gammaijUP,gcov,gcon,sqrtg)

end subroutine get_metric3plus1_both

subroutine metric3p1(x,alpha,beta,gammaijdown,gammaijUP,gcov,gcon,sqrtg)
 use metric_tools, only: get_metric
 real, intent(in)  :: x(1:3)
 real, intent(out) :: alpha,beta(1:3), gammaijdown(1:3,1:3),gammaijUP(1:3,1:3)
 real, intent(out) :: gcov(0:3,0:3),gcon(0:3,0:3),sqrtg
 real :: betaUP(1:3)
 integer :: i,j

 call get_metric(x,gcov,gcon,sqrtg)
 beta  = gcov(0,1:3)
 gammaijdown   = gcov(1:3,1:3)
 alpha = sqrt(-1./gcon(0,0))
 betaUP = gcon(0,1:3)*alpha**2
 gammaijUP = 0.
 do i=1,3
    do j=1,3
       gammaijUP(i,j) = gcon(i,j) + betaUP(i)*betaUP(j)/alpha**2
    enddo
 enddo
end subroutine metric3p1

subroutine get_u0(x,v,U0)
 use metric_tools, only: get_metric
 real, intent(in) :: x(1:3),v(1:3)
 real, intent(out) :: U0
 real :: v4(0:3), gcov(0:3,0:3), gcon(0:3,0:3), sqrtg

 call get_metric(x,gcov,gcon,sqrtg)
 v4(0) = 1.
 v4(1:3) = v(1:3)
 U0 = 1./sqrt(-dot_product_gr(v4,v4,gcov))
end subroutine get_u0

subroutine get_rderivs(position,dr)
 use metric, only: a
 real, dimension(3), intent(in) :: position
 real, dimension(3), intent(out) :: dr
 real :: x,y,z,dxdr,dydr,dzdr,r,r2sphere

 x=position(1)
 y=position(2)
 z=position(3)

 r2sphere = x**2+y**2+z**2

 r = sqrt(0.5*(r2sphere-a**2+sqrt((r2sphere-a**2))**2 + 4.*a**2*z**2))

 dxdr=x*r/(r**2+a**2)
 dydr=y*r/(r**2+a**2)
 dzdr=z/r

 dr=(/dxdr,dydr,dzdr/)

end subroutine get_rderivs

end module utils_gr
