module metric
   implicit none
   character(len=*), parameter :: metric_type = 'Kerr'

   !--- This is the Kerr metric in Kerr-Schild, 'Cartesian' coordinates

   public :: get_metric, metric_type, get_metric_derivs

   real, parameter, public :: mass1 = 1.  ! mass of central object
   real, parameter, public :: a     = 0.0 ! spin of central object
   real, parameter, public :: rs    = 2.*mass1

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
      real :: x,y,z,x2,y2,z2,r,r2,r3,r4,r_spherical,r2_spherical,a2

      a2 = a**2
      x = position(1)
      y = position(2)
      z = position(3)
      x2= x**2
      y2= y**2
      z2= z**2
      r2_spherical = x2+y2+z2
      r_spherical = sqrt(r2_spherical)
      r2 = 0.5*(r2_spherical-a**2+sqrt((r2_spherical-a**2))**2+4*a**2*z**2)
      r  = sqrt(r2)
      r3 = r*r2
      r4 = r*r3

      sqrtg = Sqrt((r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)*(a2 + r*rs)*z2)/((a2 + r2)*(r4 + a2*z2)))
      gcov(0,0) = -1 + (r3*rs)/(r4 + a2*z2)
      gcov(1,0) = (r3*rs*(r*x + a*y))/((a2 + r2)*(r4 + a2*z2))
      gcov(2,0) = (r3*rs*(-(a*x) + r*y))/((a2 + r2)*(r4 + a2*z2))
      gcov(3,0) = (r2*rs*z)/(r4 + a2*z2)
      gcov(0,1) = (r3*rs*(r*x + a*y))/((a2 + r2)*(r4 + a2*z2))
      gcov(1,1) = 1 + (r3*rs*(r*x + a*y)**2)/((a2 + r2)**2*(r4 + a2*z2))
      gcov(2,1) = (r3*rs*(r*x + a*y)*(-(a*x) + r*y))/((a2 + r2)**2*(r4 + a2*z2))
      gcov(3,1) = (r2*rs*(r*x + a*y)*z)/((a2 + r2)*(r4 + a2*z2))
      gcov(0,2) = (r3*rs*(-(a*x) + r*y))/((a2 + r2)*(r4 + a2*z2))
      gcov(1,2) = (r3*rs*(r*x + a*y)*(-(a*x) + r*y))/((a2 + r2)**2*(r4 + a2*z2))
      gcov(2,2) = 1 + (r3*rs*(a*x - r*y)**2)/((a2 + r2)**2*(r4 + a2*z2))
      gcov(3,2) = (r2*rs*(-(a*x) + r*y)*z)/((a2 + r2)*(r4 + a2*z2))
      gcov(0,3) = (r2*rs*z)/(r4 + a2*z2)
      gcov(1,3) = (r2*rs*(r*x + a*y)*z)/((a2 + r2)*(r4 + a2*z2))
      gcov(2,3) = (r2*rs*(-(a*x) + r*y)*z)/((a2 + r2)*(r4 + a2*z2))
      gcov(3,3) = (r4 + (a2 + r*rs)*z2)/(r4 + a2*z2)
      gcon(0,0) = -((r3*(a2*r + r3 + rs*(x2 + y2)) + (a2 + r2)*(a2 + r*rs)*z2)/(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2))&
      & + (a2 + r2)*(a2 + r*rs)*z2))
      gcon(1,0) = (r3*rs*(r*x + a*y))/(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)*(a2 + r*rs)*z2)
      gcon(2,0) = (r3*rs*(-(a*x) + r*y))/(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)*(a2 + r*rs)*z2)
      gcon(3,0) = (r2*(a2 + r2)*rs*z)/(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)*(a2 + r*rs)*z2)
      gcon(0,1) = (r3*rs*(r*x + a*y))/(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)*(a2 + r*rs)*z2)
      gcon(1,1) = 1 - (r3*rs*(r*x + a*y)**2)/((a2 + r2)*(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)*(a2 + r*rs)&
      &*z2))
      gcon(2,1) = -((r3*rs*(r*x + a*y)*(-(a*x) + r*y))/((a2 + r2)*(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)&
      &*(a2 + r*rs)*z2)))
      gcon(3,1) = -((r2*rs*(r*x + a*y)*z)/(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)*(a2 + r*rs)*z2))
      gcon(0,2) = (r3*rs*(-(a*x) + r*y))/(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)*(a2 + r*rs)*z2)
      gcon(1,2) = -((r3*rs*(r*x + a*y)*(-(a*x) + r*y))/((a2 + r2)*(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)&
      &*(a2 + r*rs)*z2)))
      gcon(2,2) = 1 - (r3*rs*(a*x - r*y)**2)/((a2 + r2)*(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)*(a2 + r*rs)&
      &*z2))
      gcon(3,2) = (r2*rs*(a*x - r*y)*z)/(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)*(a2 + r*rs)*z2)
      gcon(0,3) = (r2*(a2 + r2)*rs*z)/(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)*(a2 + r*rs)*z2)
      gcon(1,3) = -((r2*rs*(r*x + a*y)*z)/(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)*(a2 + r*rs)*z2))
      gcon(2,3) = (r2*rs*(a*x - r*y)*z)/(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + (a2 + r2)*(a2 + r*rs)*z2)
      gcon(3,3) = (r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2)) + a2*(a2 + r2)*z2)/(r3*(r3 + a2*(r - rs) - r2*rs + rs*(x2 + y2))&
      & + (a2 + r2)*(a2 + r*rs)*z2)
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
      real :: x,y,z,x2,y2,z2,z4,r,r2,r3,r4,r5,r6,r7,r8,a2,a3,a4,a5
      real :: r_spherical, r2_spherical

      a2 = a**2
      a3 = a*a2
      a4 = a*a3
      a5 = a*a4
      x = position(1)
      y = position(2)
      z = position(3)
      x2= x**2
      y2= y**2
      z2= z**2
      z4= z2*z2
      r2_spherical = x2+y2+z2
      r_spherical = sqrt(r2_spherical)
      r2 = 0.5*(r2_spherical-a**2+sqrt((r2_spherical-a**2))**2+4*a**2*z**2)
      r  = sqrt(r2)
      r3 = r*r2
      r4 = r*r3
      r5 = r*r4
      r6 = r*r5
      r7 = r*r6
      r8 = r*r7

      dgcovdx(0,0) = -((r5*rs*x*(r4 - 3*a2*z2))/(r4 + a2*z2)**3)
      dgcovdx(1,0) = (r4*rs*(r8*(a2 + r2 - 2*x2) - (a3 + 3*a*r2)*r5*x*y + a2*r*(2*r3*(r2 + x2) + 2*a2*(r3 + 2*r*x2)&
      & + a*(3*a2 + r2)*x*y)*z2 + a4*(a2 + r2)*z4))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdx(2,0) = (r3*rs*(-(r6*(a3*(r - x)*(r + x) + a*(r4 - 3*r2*x2) + 2*r3*x*y)) - a2*r2*(2*a*r2*(a2 + r2)&
       + a*(3*a2 + r2)*x2 - 2*r*(2*a2 + r2)*x*y)*z2 - a5*(a2 + r2)*z4))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdx(3,0) = (2*rs*x*z*(-r8 + a2*r4*z2))/(r4 + a2*z2)**3
      dgcovdx(0,1) = (r4*rs*(r8*(a2 + r2 - 2*x2) - (a3 + 3*a*r2)*r5*x*y + a2*r*(2*r3*(r2 + x2) + 2*a2*(r3 + 2*r*x2)&
      & + a*(3*a2 + r2)*x*y)*z2 + a4*(a2 + r2)*z4))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdx(1,1) = (r4*rs*(r*x + a*y)*(r5*(2*r5 - 3*r3*x2 + a2*r*(2*r2 + x2) - a3*x*y - 5*a*r2*x*y) + a2*r*(r3*(4*r2 + x2)&
      & + a2*(4*r3 + 5*r*x2) + 3*a3*x*y - a*r2*x*y)*z2 + 2*a4*(a2 + r2)*z4))/((a2 + r2)**3*(r4 + a2*z2)**3)
      dgcovdx(2,1) = (r3*rs*(4*(a2 + r2)*r6*x*(r*x + a*y)*(a*x - r*y) + r2*x*(-3*a4*x*y + 6*a2*r2*x*y + r4*x*y&
      & + 4*a3*r*(-x2 + y2))*(r4 + a2*z2) + (a2 + r2)*(-2*a*r*x - a2*y + r2*y)*(r4 + a2*z2)**2))/((a2 + r2)**3*(r4 + a2*z2)**3)
      dgcovdx(3,1) = (r3*rs*z*(-4*(a2 + r2)*r5*x*(r*x + a*y) + r*x*(3*a2*r*x + r3*x + 2*a3*y)*(r4 + a2*z2)&
      & + (a2 + r2)*(r4 + a2*z2)**2))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdx(0,2) = (r3*rs*(-(r6*(a3*(r - x)*(r + x) + a*(r4 - 3*r2*x2) + 2*r3*x*y)) - a2*r2*(2*a*r2*(a2 + r2)&
      & + a*(3*a2 + r2)*x2 - 2*r*(2*a2 + r2)*x*y)*z2 - a5*(a2 + r2)*z4))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdx(1,2) = (r3*rs*(4*(a2 + r2)*r6*x*(r*x + a*y)*(a*x - r*y) + r2*x*(-3*a4*x*y + 6*a2*r2*x*y + r4*x*y&
      & + 4*a3*r*(-x2 + y2))*(r4 + a2*z2) + (a2 + r2)*(-2*a*r*x - a2*y + r2*y)*(r4 + a2*z2)**2))/((a2 + r2)**3*(r4 + a2*z2)**3)
      dgcovdx(2,2) = (r3*rs*(a*x - r*y)*(r6*(a3*(2*r2 - x2) + a*(2*r4 - 5*r2*x2) - a2*r*x*y + 3*r3*x*y)&
      & + a2*r2*(4*a*r2*(a2 + r2) + a*(3*a2 - r2)*x2 - r*(5*a2 + r2)*x*y)*z2 + 2*a5*(a2 + r2)*z4))/((a2 + r2)**3*(r4 + a2*z2)**3)
      dgcovdx(3,2) = (r2*rs*z*(-(r6*(a*r2*(a2 + r2) - 2*a*(a2 + 2*r2)*x2 + r*(a2 + 3*r2)*x*y)) + a2*r2*(-2*a*(r4 + a2*(r2 + x2))&
      & + r*(3*a2 + r2)*x*y)*z2 - a5*(a2 + r2)*z4))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdx(0,3) = (2*rs*x*z*(-r8 + a2*r4*z2))/(r4 + a2*z2)**3
      dgcovdx(1,3) = (r3*rs*z*(-4*(a2 + r2)*r5*x*(r*x + a*y) + r*x*(3*a2*r*x + r3*x + 2*a3*y)*(r4 + a2*z2)&
      & + (a2 + r2)*(r4 + a2*z2)**2))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdx(2,3) = (r2*rs*z*(-(r6*(a*r2*(a2 + r2) - 2*a*(a2 + 2*r2)*x2 + r*(a2 + 3*r2)*x*y)) + a2*r2*(-2*a*(r4 + a2*(r2 + x2))&
      & + r*(3*a2 + r2)*x*y)*z2 - a5*(a2 + r2)*z4))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdx(3,3) = (r3*rs*x*z2*(-3*r4 + a2*z2))/(r4 + a2*z2)**3
      dgcovdy(0,0) = -((r5*rs*y*(r4 - 3*a2*z2))/(r4 + a2*z2)**3)
      dgcovdy(1,0) = (r3*rs*(-2*r**9*x*y + a3*r6*(r - y)*(r + y) + a*r8*(r2 - 3*y2) + a2*r2*(2*a*r2*(a2 + r2)&
      & + 2*r*(2*a2 + r2)*x*y + a*(3*a2 + r2)*y2)*z2 + a5*(a2 + r2)*z4))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdy(2,0) = (r4*rs*((a2 + r2)*r8 + (a3 + 3*a*r2)*r5*x*y - 2*r8*y2 + a2*r*(2*(a2 + r2)*r3 - a*(3*a2 + r2)*x*y&
      & + 2*r*(2*a2 + r2)*y2)*z2 + a4*(a2 + r2)*z4))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdy(3,0) = (2*rs*y*z*(-r8 + a2*r4*z2))/(r4 + a2*z2)**3
      dgcovdy(0,1) = (r3*rs*(-2*r**9*x*y + a3*r6*(r - y)*(r + y) + a*r8*(r2 - 3*y2) + a2*r2*(2*a*r2*(a2 + r2)&
      & + 2*r*(2*a2 + r2)*x*y + a*(3*a2 + r2)*y2)*z2 + a5*(a2 + r2)*z4))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdy(1,1) = (r3*rs*(r*x + a*y)*(r6*(2*a*r2*(a2 + r2) + r*(a2 - 3*r2)*x*y - a*(a2 + 5*r2)*y2)&
      & + a2*r2*(4*a*r2*(a2 + r2) + r*(5*a2 + r2)*x*y + a*(3*a2 - r2)*y2)*z2 + 2*a5*(a2 + r2)*z4))&
      &/((a2 + r2)**3*(r4 + a2*z2)**3)
      dgcovdy(2,1) = (r3*rs*(4*(a2 + r2)*r6*y*(r*x + a*y)*(a*x - r*y) + r2*y*(-3*a4*x*y + 6*a2*r2*x*y&
      & + r4*x*y + 4*a3*r*(-x2 + y2))*(r4 + a2*z2) + (a2 + r2)*(-(a2*x) + r2*x + 2*a*r*y)*(r4 + a2*z2)**2))&
      &/((a2 + r2)**3*(r4 + a2*z2)**3)
      dgcovdy(3,1) = (rs*z*(-4*(a2 + r2)*r8*y*(r*x + a*y) + r4*y*(3*a2*r*x + r3*x + 2*a3*y)*(r4 + a2*z2)&
      & + a*(a2 + r2)*(r5 + a2*r*z2)**2))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdy(0,2) = (r4*rs*((a2 + r2)*r8 + (a3 + 3*a*r2)*r5*x*y - 2*r8*y2 + a2*r*(2*(a2 + r2)*r3 - a*(3*a2 + r2)*x*y&
      & + 2*r*(2*a2 + r2)*y2)*z2 + a4*(a2 + r2)*z4))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdy(1,2) = (r3*rs*(4*(a2 + r2)*r6*y*(r*x + a*y)*(a*x - r*y) + r2*y*(-3*a4*x*y + 6*a2*r2*x*y&
      & + r4*x*y + 4*a3*r*(-x2 + y2))*(r4 + a2*z2) + (a2 + r2)*(-(a2*x) + r2*x + 2*a*r*y)*(r4 + a2*z2)**2))&
      &/((a2 + r2)**3*(r4 + a2*z2)**3)
      dgcovdy(2,2) = (r4*rs*(-(a*x) + r*y)*(r5*(2*(a2 + r2)*r3 + a*(a2 + 5*r2)*x*y + r*(a2 - 3*r2)*y2)&
      & + a2*r*(-3*a3*x*y + a*r2*x*y + r3*(4*r2 + y2) + a2*(4*r3 + 5*r*y2))*z2 + 2*a4*(a2 + r2)*z4))&
      &/((a2 + r2)**3*(r4 + a2*z2)**3)
      dgcovdy(3,2) = (r3*rs*z*(r5*(r5 + 2*a3*x*y + 4*a*r2*x*y + a2*r*(r - y)*(r + y) - 3*r3*y2)&
      & + a2*r*(2*(a2 + r2)*r3 - 2*a3*x*y + r*(3*a2 + r2)*y2)*z2 + a4*(a2 + r2)*z4))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdy(0,3) = (2*rs*y*z*(-r8 + a2*r4*z2))/(r4 + a2*z2)**3
      dgcovdy(1,3) = (rs*z*(-4*(a2 + r2)*r8*y*(r*x + a*y) + r4*y*(3*a2*r*x + r3*x + 2*a3*y)*(r4 + a2*z2)&
      & + a*(a2 + r2)*(r5 + a2*r*z2)**2))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdy(2,3) = (r3*rs*z*(r5*(r5 + 2*a3*x*y + 4*a*r2*x*y + a2*r*(r - y)*(r + y) - 3*r3*y2)&
      & + a2*r*(2*(a2 + r2)*r3 - 2*a3*x*y + r*(3*a2 + r2)*y2)*z2 + a4*(a2 + r2)*z4))/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdy(3,3) = (r3*rs*y*z2*(-3*r4 + a2*z2))/(r4 + a2*z2)**3
      dgcovdz(0,0) = (r3*rs*z*(-r6 + a4*z2 + 3*a2*r2*(-r2 + z2)))/(r4 + a2*z2)**3
      dgcovdz(1,0) = (-(r7*rs*(2*r*x + 3*a*y)*z) + a2*r3*rs*(2*r*x + a*y)*z**3)/(r4 + a2*z2)**3
      dgcovdz(2,0) = (r7*rs*(3*a*x - 2*r*y)*z + a2*r3*rs*(-(a*x) + 2*r*y)*z**3)/(r4 + a2*z2)**3
      dgcovdz(3,0) = (r2*rs*(r8 - 2*(a2 + r2)*r4*z2 + a2*(a2 + 2*r2)*z4))/(r4 + a2*z2)**3
      dgcovdz(0,1) = (-(r7*rs*(2*r*x + 3*a*y)*z) + a2*r3*rs*(2*r*x + a*y)*z**3)/(r4 + a2*z2)**3
      dgcovdz(1,1) = -((r3*rs*(r*x + a*y)*z*(r4*(a2*r*x + 3*r3*x + 3*a3*y + 5*a*r2*y)&
      & - a2*(3*a2*r*x + r3*x + a3*y - a*r2*y)*z2))/((a2 + r2)**2*(r4 + a2*z2)**3))
      dgcovdz(2,1) = (r7*rs*(3*a4*x*y + 4*a2*r2*x*y - 3*r4*x*y + 2*a3*r*(x - y)*(x + y)&
      & + 4*a*r3*(x - y)*(x + y))*z + a2*r3*rs*(-(a4*x*y) + 4*a2*r2*x*y + r4*x*y + 2*a3*r*(-x2 + y2))*z**3)&
      &/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdz(3,1) = (r2*rs*(r8*(r*x + a*y) - r4*(a2*r*x + 3*r3*x + 2*a3*y + 4*a*r2*y)*z2 + a2*(2*a2*r*x + r3*x + a3*y)*z4))&
      &/((a2 + r2)*(r4 + a2*z2)**3)
      dgcovdz(0,2) = (r7*rs*(3*a*x - 2*r*y)*z + a2*r3*rs*(-(a*x) + 2*r*y)*z**3)/(r4 + a2*z2)**3
      dgcovdz(1,2) = (r7*rs*(3*a4*x*y + 4*a2*r2*x*y - 3*r4*x*y + 2*a3*r*(x - y)*(x + y) + 4*a*r3*(x - y)*(x + y))*z&
      & + a2*r3*rs*(-(a4*x*y) + 4*a2*r2*x*y + r4*x*y + 2*a3*r*(-x2 + y2))*z**3)/((a2 + r2)**2*(r4 + a2*z2)**3)
      dgcovdz(2,2) = -((r3*rs*(-(a*x) + r*y)*z*(r4*(-3*a3*x - 5*a*r2*x + a2*r*y + 3*r3*y) + a2*(a*(a - r)*(a + r)*x&
      & - r*(3*a2 + r2)*y)*z2))/((a2 + r2)**2*(r4 + a2*z2)**3))
      dgcovdz(3,2) = (r2*rs*(r8*(-(a*x) + r*y) + r4*(2*a*(a2 + 2*r2)*x - r*(a2 + 3*r2)*y)*z2&
      & + a2*(-(a3*x) + 2*a2*r*y + r3*y)*z4))/((a2 + r2)*(r4 + a2*z2)**3)
      dgcovdz(0,3) = (r2*rs*(r8 - 2*(a2 + r2)*r4*z2 + a2*(a2 + 2*r2)*z4))/(r4 + a2*z2)**3
      dgcovdz(1,3) = (r2*rs*(r8*(r*x + a*y) - r4*(a2*r*x + 3*r3*x + 2*a3*y + 4*a*r2*y)*z2&
      & + a2*(2*a2*r*x + r3*x + a3*y)*z4))/((a2 + r2)*(r4 + a2*z2)**3)
      dgcovdz(2,3) = (r2*rs*(r8*(-(a*x) + r*y) + r4*(2*a*(a2 + 2*r2)*x - r*(a2 + 3*r2)*y)*z2&
      & + a2*(-(a3*x) + 2*a2*r*y + r3*y)*z4))/((a2 + r2)*(r4 + a2*z2)**3)
      dgcovdz(3,3) = (rs*z*(2*r**9 - (a2 + 3*r2)*r5*z2 + a2*r*(a2 + r2)*z4))/(r4 + a2*z2)**3


   end subroutine get_metric_derivs

end module metric
