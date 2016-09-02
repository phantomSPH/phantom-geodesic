module metric
 implicit none
 character(len=*), parameter :: metric_type = 'Schwarzschild'

 public :: get_metric, metric_type, get_metric_derivs

 real, parameter, public :: mass1 = 1.  ! mass of central object

 private

contains

 !----------------------------------------------------------------
 !+
 !  Compute the metric tensor in both covariant (gcov) and
 !  contravariant (gcon) form
 !+
 !----------------------------------------------------------------
 pure subroutine get_metric(position,gcov,gcon,sqrtg)
  real,    intent(in)  :: position(3)
  real,    intent(out) :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
  real :: rs,r,r2,r3,rs_on_r3,coeff,x,y,z,x2,y2,z2,term

  gcov = 0.
  gcon = 0.
  rs = 2.*mass1
  r2 = dot_product(position,position)
  r  = sqrt(r2)
  r3 = r*r2
  rs_on_r3 = rs/r3
  x  = position(1)
  y  = position(2)
  z  = position(3)
  x2 = x**2
  y2 = y**2
  z2 = z**2

   term  = (1.-rs/r)
   coeff = 1./term
   gcov(0,0) = -term
   gcov(1,1) = coeff*(1.-rs_on_r3*(y2+z2))
   gcov(2,2) = coeff*(1.-rs_on_r3*(x2+z2))
   gcov(3,3) = coeff*(1.-rs_on_r3*(x2+y2))
   gcov(1,2) = coeff*x*y*rs_on_r3
   gcov(2,1) = gcov(1,2)
   gcov(1,3) = coeff*x*z*rs_on_r3
   gcov(3,1) = gcov(1,3)
   gcov(2,3) = coeff*y*z*rs_on_r3
   gcov(3,2) = gcov(2,3)

   sqrtg=1.

   gcon(0,0) = -1./term
   gcon(1,1) = 1.-rs_on_r3*x2
   gcon(2,2) = 1.-rs_on_r3*y2
   gcon(3,3) = 1.-rs_on_r3*z2
   gcon(1,2) = -rs_on_r3*x*y
   gcon(2,1) = gcon(1,2)
   gcon(1,3) = -rs_on_r3*x*z
   gcon(3,1) = gcon(1,3)
   gcon(2,3) = -rs_on_r3*y*z
   gcon(3,2) = gcon(2,3)

 end subroutine get_metric


 !----------------------------------------------------------------
 !+
 !  Compute the derivatives of the metric tensor in both covariant
 !  (dgcovdx, dgcovdy, dgcovdz) and contravariant (dgcondx, dgcondy,
 !  dgcondz) form
 !+
 !----------------------------------------------------------------
 pure subroutine get_metric_derivs(position,dgcovdx, dgcovdy, dgcovdz)
  real,    intent(in)  :: position(3)
  real,    intent(out) :: dgcovdx(0:3,0:3), dgcovdy(0:3,0:3), dgcovdz(0:3,0:3)
  ! real,    intent(out) :: dgcondx(0:3,0:3), dgcondy(0:3,0:3), dgcondz(0:3,0:3)
  real :: x,y,z, rs,r,r2,r3,r4,r5,rs_on_r3,x2,y2,z2,rs2

  dgcovdx = 0.
  dgcovdy = 0.
  dgcovdz = 0.
  ! dgcondx = 0.
  ! dgcondy = 0.
  ! dgcondz = 0.

  x = position(1)
  y = position(2)
  z = position(3)
  x2= x**2
  y2= y**2
  z2= z**2
  rs = 2.*mass1
  r2 = dot_product(position,position)
  r  = sqrt(r2)
  r3 = r*r2
  r4 = r2*r2
  r5 = r*r4
  rs_on_r3 = rs/r3
  rs2 = rs**2

   ! dx
   dgcovdx(0,0) = -((rs*x)/r3)
   dgcovdx(0,1) = 0
   dgcovdx(0,2) = 0
   dgcovdx(0,3) = 0
   dgcovdx(1,0) = 0
   dgcovdx(1,1) = (x*(2*r4 - 2*r3*rs - 2*rs2*(y2 + z2) - 2*r4&
    + r*rs*(x2 + 4*(y2 + z2))))/(r4*(r - rs)**2)
   dgcovdx(1,2) = (rs*(r3 - r2*rs - 3*r*x2 + 2*rs*x2)*y)/(r4*(r - rs)**2)
   dgcovdx(1,3) = (rs*(r3 - r2*rs - 3*r*x2 + 2*rs*x2)*z)/(r4*(r - rs)**2)
   dgcovdx(2,0) = 0
   dgcovdx(2,1) = (rs*(r3 - r2*rs - 3*r*x2 + 2*rs*x2)*y)/(r4*(r - rs)**2)
   dgcovdx(2,2) = (x*(2*r4 - 4*r3*rs + 2*r2*(rs2 - r2) &
   - 2*rs2*(x2 + z2) + r*rs*(4*x2 + y2 + 4*z2)))/(r4*(r - rs)**2)
   dgcovdx(2,3) = (rs*(-3*r + 2*rs)*x*y*z)/(r4*(r - rs)**2)
   dgcovdx(3,0) = 0
   dgcovdx(3,1) = (rs*(r3 - r2*rs - 3*r*x2 + 2*rs*x2)*z)/(r4*(r - rs)**2)
   dgcovdx(3,2) = (rs*(-3*r + 2*rs)*x*y*z)/(r4*(r - rs)**2)
   dgcovdx(3,3) = (-2*(r - rs)**2*x*(-z2) + r*(-2*r + rs)*x*z2)&
   /(r4*(r - rs)**2)

   ! dgcondx(0,0) = (rs*x)/(r*(r - rs)**2)
   ! dgcondx(0,1) = 0
   ! dgcondx(0,2) = 0
   ! dgcondx(0,3) = 0
   ! dgcondx(1,0) = 0
   ! dgcondx(1,1) = -((x*(-2*r3 + 2*r2*rs - 3*rs*x2 + 2*r3))/r5)
   ! dgcondx(1,2) = -((rs*(r2 - 3*x2)*y)/r5)
   ! dgcondx(1,3) = -((rs*(r2 - 3*x2)*z)/r5)
   ! dgcondx(2,0) = 0
   ! dgcondx(2,1) = -((rs*(r2 - 3*x2)*y)/r5)
   ! dgcondx(2,2) = (x*(2*r3 + 3*rs*y2 - 2*r3))/r5
   ! dgcondx(2,3) = (3*rs*x*y*z)/r5
   ! dgcondx(3,0) = 0
   ! dgcondx(3,1) = -((rs*(r2 - 3*x2)*z)/r5)
   ! dgcondx(3,2) = (3*rs*x*y*z)/r5
   ! dgcondx(3,3) = (x*(2*r3 + 3*rs*z2 - 2*r3))/r5

   ! dy
   dgcovdy(0,0) = -((rs*y)/r3)
   dgcovdy(0,1) = 0
   dgcovdy(0,2) = 0
   dgcovdy(0,3) = 0
   dgcovdy(1,0) = 0
   dgcovdy(1,1) = (y*(2*r4 - 4*r3*rs + 2*r2*(rs2 - r2) &
   - 2*rs2*(y2 + z2) + r*rs*(x2 + 4*(y2 + z2))))/(r4*(r - rs)**2)
   dgcovdy(1,2) = (rs*x*(r3 - r2*rs - 3*r*y2 + 2*rs*y2))/(r4*(r - rs)**2)
   dgcovdy(1,3) = (rs*(-3*r + 2*rs)*x*y*z)/(r4*(r - rs)**2)
   dgcovdy(2,0) = 0
   dgcovdy(2,1) = (rs*x*(r3 - r2*rs - 3*r*y2 + 2*rs*y2))/(r4*(r - rs)**2)
   dgcovdy(2,2) = (y*(2*r4 - 2*r3*rs - 2*rs2*(x2 + z2) - 2*r4&
    + r*rs*(4*x2 + y2 + 4*z2)))/(r4*(r - rs)**2)
   dgcovdy(2,3) = (rs*(r3 - r2*rs - 3*r*y2 + 2*rs*y2)*z)/(r4*(r - rs)**2)
   dgcovdy(3,0) = 0
   dgcovdy(3,1) = (rs*(-3*r + 2*rs)*x*y*z)/(r4*(r - rs)**2)
   dgcovdy(3,2) = (rs*(r3 - r2*rs - 3*r*y2 + 2*rs*y2)*z)/(r4*(r - rs)**2)
   dgcovdy(3,3) = (-2*(r - rs)**2*y*(-z2) + r*(-2*r + rs)*y*z2)&
   /(r4*(r - rs)**2)

   ! dgcondy(0,0) = (rs*y)/(r*(r - rs)**2)
   ! dgcondy(0,1) = 0
   ! dgcondy(0,2) = 0
   ! dgcondy(0,3) = 0
   ! dgcondy(1,0) = 0
   ! dgcondy(1,1) = (y*(2*r3 + 3*rs*x2 - 2*r3))/r5
   ! dgcondy(1,2) = -((rs*x*(r2 - 3*y2))/r5)
   ! dgcondy(1,3) = (3*rs*x*y*z)/r5
   ! dgcondy(2,0) = 0
   ! dgcondy(2,1) = -((rs*x*(r2 - 3*y2))/r5)
   ! dgcondy(2,2) = -((y*(-2*r3 + 2*r2*rs - 3*rs*y2 + 2*r3))/r5)
   ! dgcondy(2,3) = -((rs*(r2 - 3*y2)*z)/r5)
   ! dgcondy(3,0) = 0
   ! dgcondy(3,1) = (3*rs*x*y*z)/r5
   ! dgcondy(3,2) = -((rs*(r2 - 3*y2)*z)/r5)
   ! dgcondy(3,3) = (y*(2*r3 + 3*rs*z2 - 2*r3))/r5

   ! dz
   dgcovdz(0,0) = -((rs*z)/r3)
   dgcovdz(0,1) = 0
   dgcovdz(0,2) = 0
   dgcovdz(0,3) = 0
   dgcovdz(1,0) = 0
   dgcovdz(1,1) = (z*(2*r4 - 4*r3*rs + 2*r2*(rs2 - r2) &
   - 2*rs2*(y2 + z2) + r*rs*(x2 + 4*(y2 + z2))))/(r4*(r - rs)**2)
   dgcovdz(1,2) = (rs*(-3*r + 2*rs)*x*y*z)/(r4*(r - rs)**2)
   dgcovdz(1,3) = (rs*x*(r3 - r2*rs - 3*r*z2 + 2*rs*z2))/(r4*(r - rs)**2)
   dgcovdz(2,0) = 0
   dgcovdz(2,1) = (rs*(-3*r + 2*rs)*x*y*z)/(r4*(r - rs)**2)
   dgcovdz(2,2) = (z*(2*r4 - 4*r3*rs + 2*r2*(rs2 - r2) -&
   2*rs2*(x2 + z2) + r*rs*(4*x2 + y2 + 4*z2)))/(r4*(r - rs)**2)
   dgcovdz(2,3) = (rs*y*(r3 - r2*rs - 3*r*z2 + 2*rs*z2))/(r4*(r - rs)**2)
   dgcovdz(3,0) = 0
   dgcovdz(3,1) = (rs*x*(r3 - r2*rs - 3*r*z2 + 2*rs*z2))/(r4*(r - rs)**2)
   dgcovdz(3,2) = (rs*y*(r3 - r2*rs - 3*r*z2 + 2*rs*z2))/(r4*(r - rs)**2)
   dgcovdz(3,3) = (z*(2*(r - rs)*(r3 - r*(x2 + y2) + rs*(x2 + y2)) + &
   r*(-2*r + rs)*z2))/(r4*(r - rs)**2)

   ! dgcondz(0,0) = (rs*z)/(r*(r - rs)**2)
   ! dgcondz(0,1) = 0
   ! dgcondz(0,2) = 0
   ! dgcondz(0,3) = 0
   ! dgcondz(1,0) = 0
   ! dgcondz(1,1) = (z*(2*r3 + 3*rs*x2 - 2*r3))/r5
   ! dgcondz(1,2) = (3*rs*x*y*z)/r5
   ! dgcondz(1,3) = -((rs*x*(r2 - 3*z2))/r5)
   ! dgcondz(2,0) = 0
   ! dgcondz(2,1) = (3*rs*x*y*z)/r5
   ! dgcondz(2,2) = (z*(2*r3 + 3*rs*y2 - 2*r3))/r5
   ! dgcondz(2,3) = -((rs*y*(r2 - 3*z2))/r5)
   ! dgcondz(3,0) = 0
   ! dgcondz(3,1) = -((rs*x*(r2 - 3*z2))/r5)
   ! dgcondz(3,2) = -((rs*y*(r2 - 3*z2))/r5)
   ! dgcondz(3,3) = (z*(2*r3 - 2*r2*rs + 3*rs*z2 - 2*r3))/r5

 end subroutine get_metric_derivs

end module metric
