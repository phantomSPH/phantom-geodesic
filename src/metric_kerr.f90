module metric
 implicit none
 character(len=*), parameter :: metric_type = 'Kerr'

 public :: get_metric, metric_type, get_metric_derivs

 real, parameter, public :: mass1 = 1.  ! mass of central object
 real, parameter, public :: alpha = 0.0 ! spin of central object

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
  real :: x,y,z, rs,r,r2,r3,r4,r5,r6,x2,y2,z2
  real :: alpha2,costheta,costheta2,rho2,delta,A,Adeltar6,alphars,alpha2plusr2
  real :: rrho2, Ar4, Ar4plusdeltaz2,alpha2rs,r3rho2,r2minusdelta,r2rho2,rrs

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
  r6 = r*r5
  alpha2 = alpha**2

  costheta  = z/r
  costheta2 = costheta2**2

  rho2 = r2 + alpha2*costheta2
  delta = r2 - rs*r + alpha2
  A = 1.-costheta2

  Adeltar6 = A*delta*r6
  alphars = alpha*rs
  alpha2plusr2 = alpha2 + r2
  rrho2 = r*rho2
  Ar4 = A*r4
  Ar4plusdeltaz2= Ar4 + delta*z2
  alpha2rs = alpha2*rs
  r3rho2 = r3*rho2
  r2minusdelta = r2 - delta
  r2rho2 = r2*rho2
  rrs = r*rs

  sqrtg = rho2/r2

  gcov(0,0) = -1. + (rrs)/rho2
  gcov(1,0) = (alphars*y)/(rrho2)
  gcov(2,0) = -((alphars*x)/(rrho2))
  gcov(3,0) = 0.

  gcov(0,1) = gcov(1,0)
  gcov(1,1) = (alpha2plusr2*y2)/(Ar4) + (alpha2rs*y2)/(r3rho2) + (rho2*x2*(Ar4plusdeltaz2))/(Adeltar6)
  gcov(2,1) = -((alpha2plusr2*x*y)/(Ar4)) - (alpha2rs*x*y)/(r3rho2) + (rho2*x*y*(Ar4plusdeltaz2))/(Adeltar6)
  gcov(3,1) = (r2minusdelta*rho2*x*z)/(delta*r4)

  gcov(0,2) = gcov(2,0)
  gcov(1,2) = gcov(2,1)
  gcov(2,2) = (alpha2plusr2*x2)/(Ar4) + (alpha2rs*x2)/(r3rho2) + (rho2*y2*(Ar4plusdeltaz2))/(Adeltar6)
  gcov(3,2) = (r2minusdelta*rho2*y*z)/(delta*r4)

  gcov(0,3) = 0.
  gcov(1,3) = gcov(3,1)
  gcov(2,3) = gcov(3,2)
  gcov(3,3) = ((A*delta + z2)*rho2)/(delta*r2)

  !----------

  gcon(0,0) = -1. - (r*alpha2plusr2*rs)/(delta*rho2)
  gcon(1,0) = (alpha*rrs*y)/(delta*rho2)
  gcon(2,0) = -((alpha*rrs*x)/(delta*rho2))
  gcon(3,0) = 0.

  gcon(0,1) = gcon(1,0)
  gcon(1,1) = (((rho2 - rrs)*y**2)/delta + (x**2*(A*delta + z2))/r2)/(A*rho2)
  gcon(2,1) = x*y*((alpha2*rrs)/(delta*alpha2plusr2*rho2) + (-(1/alpha2plusr2) + (A*delta + z2)/r2rho2)/A)
  gcon(3,1) = -((r2minusdelta*x*z)/r2rho2)

  gcon(0,2) = gcon(2,0)
  gcon(1,2) = gcon(2,1)
  gcon(2,2) = (((rho2 - rrs)*x**2)/delta + (y**2*(A*delta + z2))/r2)/(A*rho2)
  gcon(3,2) = -((r2minusdelta*y*z)/r2rho2)

  gcon(0,3) = 0.
  gcon(1,3) = gcon(3,1)
  gcon(2,3) = gcon(3,2)
  gcon(3,3) = (r4 - r2minusdelta*z2)/r2rho2

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
    
 end subroutine get_metric_derivs

end module metric
