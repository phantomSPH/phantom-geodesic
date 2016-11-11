module metric
 implicit none
 character(len=*), parameter :: metric_type = 'Minkowski'

 real, parameter, public :: mass1 = 0.  ! mass of central object
 real, parameter, public :: a     = 0.  ! spin of central object
 real, parameter, public :: rs    = 2.*mass1

contains

 !----------------------------------------------------------------
 !+
 !  Compute the metric tensor in both covariant (gcov) and
 !  contravariant (gcon) form
 !+
 !----------------------------------------------------------------
 subroutine get_metric_cartesian(position,gcov,gcon,sqrtg)
  real,    intent(in)  :: position(3)
  real,    intent(out) :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
  gcov = 0.
  gcon = 0.
  if (metric_type .eq. 'Minkowski') then
   gcov(0,0) = -1.
   gcov(1,1) = 1.
   gcov(2,2) = 1.
   gcov(3,3) = 1.
   sqrtg=1.
   gcon = gcov
  endif
 end subroutine get_metric_cartesian

 subroutine get_metric_spherical(position,gcov,gcon,sqrtg)
    real,    intent(in)  :: position(3)
    real,    intent(out) :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
    STOP 'metric in spherical coordinates not implemented'
 end subroutine get_metric_spherical

 subroutine metric_cartesian_derivatives(position,dgcovdx, dgcovdy, dgcovdz)
    real,    intent(in)  :: position(3)
    real,    intent(out) :: dgcovdx(0:3,0:3), dgcovdy(0:3,0:3), dgcovdz(0:3,0:3)
    dgcovdx = 0.
    dgcovdy = 0.
    dgcovdz = 0.
end subroutine metric_cartesian_derivatives

 subroutine metric_spherical_derivatives(position,dgcovdr, dgcovdtheta, dgcovdphi)
   real, intent(in) :: position(3)
   real, intent(out), dimension(0:3,0:3) :: dgcovdr,dgcovdtheta,dgcovdphi
   real :: r, theta
      STOP 'no sphericals metric derivs implemented'

 end subroutine metric_spherical_derivatives

 subroutine cartesian2spherical(xcart,xspher)
  real, intent(in) :: xcart(3)
  real, intent(out) ::xspher(3)
  real :: x,y,z
  real :: r,theta,phi
      stop "no cartesian2spherical implemented"
 end subroutine cartesian2spherical

 subroutine get_jacobian(position,dxdx)
  real, intent(in), dimension(3) :: position
  real, intent(out), dimension(0:3,0:3) :: dxdx
  stop "no Jacobian implemented"
 end subroutine get_jacobian
end module metric
