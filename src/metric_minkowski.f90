module metric
 implicit none
 character(len=*), parameter :: metric_type = 'Minkowski'

 real, parameter, public :: mass1 = 1000.  ! mass of central object
 real,            public :: a     = 0.  ! spin of central object
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

 gcov = 0.*position(1)
 gcon = 0.

 gcov(0,0) = -1.
 gcov(1,1) = 1.
 gcov(2,2) = 1.
 gcov(3,3) = 1.

 gcon      = gcov

 sqrtg     = 1.
end subroutine get_metric_cartesian

subroutine get_metric_spherical(position,gcov,gcon,sqrtg)
 real,    intent(in)  :: position(3)
 real,    intent(out) :: gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
 real :: r2,sintheta

 gcov = 0.
 gcon = 0.

 r2       = position(1)**2
 sintheta = sin(position(2))

 gcov(0,0) = -1.
 gcov(1,1) = 1.
 gcov(2,2) = r2
 gcov(3,3) = r2*sintheta**2

 gcon(0,0) = -1.
 gcon(1,1) = 1.
 gcon(2,2) = 1./r2
 gcov(3,3) = 1./gcov(3,3)

 sqrtg     = r2*sintheta
end subroutine get_metric_spherical

subroutine metric_cartesian_derivatives(position,dgcovdx, dgcovdy, dgcovdz)
 real,    intent(in)  :: position(3)
 real,    intent(out) :: dgcovdx(0:3,0:3), dgcovdy(0:3,0:3), dgcovdz(0:3,0:3)
 dgcovdx = 0.*position(1)
 dgcovdy = 0.
 dgcovdz = 0.
end subroutine metric_cartesian_derivatives

subroutine metric_spherical_derivatives(position,dgcovdr, dgcovdtheta, dgcovdphi)
 real, intent(in) :: position(3)
 real, intent(out), dimension(0:3,0:3) :: dgcovdr,dgcovdtheta,dgcovdphi
 real :: r, theta

 r     = position(1)
 theta = position(2)

 dgcovdr     = 0.
 dgcovdtheta = 0.
 dgcovdphi   = 0.

 dgcovdr(2,2) = 2*r
 dgcovdr(3,3) = 2*r*sin(theta)**2

 dgcovdtheta(3,3) = 2*r**2*cos(theta)*sin(theta)

end subroutine metric_spherical_derivatives

subroutine cartesian2spherical(xcart,xspher,vcart,vspher)
 real, intent(in)  :: xcart(3)
 real, intent(out) :: xspher(3)
 real, intent(in),  optional :: vcart(3)
 real, intent(out), optional :: vspher(3)
 real :: x,y,z
 real :: r,theta,phi
 real :: vx,vy,vz,vr,x2,y2

 if (present(vcart).and..not.present(vspher)) print*,'Warning: cartesian2spherical not converting velocity'
 if (present(vspher).and..not.present(vcart)) print*,'Warning: cartesian2spherical not converting velocity'

 x  = xcart(1)
 y  = xcart(2)
 z  = xcart(3)
 x2 = x**2
 y2 = y**2

 r  = sqrt(x2+y2+z**2)
 theta = acos(z/r)
 phi   = atan2(y,x)

 xspher   = (/r,theta,phi/)

 if (present(vcart).and.present(vspher)) then
   vx        = vcart(1)
   vy        = vcart(2)
   vz        = vcart(3)
   vr        = (x*vx + y*vy + z*vz)/r
   vspher(1) = vr
   vspher(2) = (z*vr - vz*r)/(r*sqrt(x2+y2))
   vspher(3) = (x*vy - y*vx)/(x2+y2)
 endif

end subroutine cartesian2spherical

subroutine spherical2cartesian(xspher,xcart,vspher,vcart)
 real, intent(in)  :: xspher(3)
 real, intent(out) :: xcart(3)
 real, intent(in),  optional :: vspher(3)
 real, intent(out), optional :: vcart(3)
 real :: x,y,z,r,theta,phi
 real :: sintheta,costheta,sinphi,cosphi,vr,vtheta,vphi

 if (present(vcart).and..not.present(vspher)) print*,'Warning: spherical2cartesian not converting velocity'
 if (present(vspher).and..not.present(vcart)) print*,'Warning: spherical2cartesian not converting velocity'

 r     = xspher(1)
 theta = xspher(2)
 phi   = xspher(3)

 sintheta = sin(theta)
 costheta = cos(theta)
 sinphi   = sin(phi)
 cosphi   = cos(phi)

 x = r*sintheta*cosphi
 y = r*sintheta*sinphi
 z = r*costheta

 xcart = (/x,y,z/)

 if (present(vspher).and.present(vcart)) then
    vr     = vspher(1)
    vtheta = vspher(2)
    vphi   = vspher(3)
    vcart(1) = sintheta*cosphi*vr + r*(costheta*cosphi*vtheta - sintheta*sinphi*vphi)
    vcart(2) = sintheta*sinphi*vr + r*(costheta*sinphi*vtheta + sintheta*cosphi*vphi)
    vcart(3) = costheta*vr - r*sintheta*vtheta
 endif

end subroutine spherical2cartesian

subroutine get_jacobian(position,dxdx)
 real, intent(in), dimension(3) :: position
 real, intent(out), dimension(0:3,0:3) :: dxdx
 real, dimension(3) :: dSPHERICALdx,dSPHERICALdy,dSPHERICALdz
 real :: drdx,drdy,drdz
 real :: dthetadx,dthetady,dthetadz
 real :: dphidx,dphidy,dphidz
 real :: x,y,z,x2,y2,z2,r2,r,rcyl2,rcyl
 x  = position(1)
 y  = position(2)
 z  = position(3)
 x2 = x**2
 y2 = y**2
 z2 = z**2
 r2 = x2+y2+z2
 r  = sqrt(r2)
 rcyl2 = x2+y2
 rcyl  = sqrt(x2+y2)

 drdx = x/r
 drdy = y/r
 drdz = z/r

 dthetadx = x*z/(r2*rcyl)
 dthetady = y*z/(r2*rcyl)
 dthetadz = -rcyl/r2

 dphidx = -y/(x2+y2)
 dphidy = x/(x2+y2)
 dphidz = 0.

 dSPHERICALdx=(/drdx,dthetadx,dphidx/)
 dSPHERICALdy=(/drdy,dthetady,dphidy/)
 dSPHERICALdz=(/drdz,dthetadz,dphidz/)

 dxdx        = 0.
 dxdx(0,0)   = 1.
 dxdx(1:3,1) = dSPHERICALdx
 dxdx(1:3,2) = dSPHERICALdy
 dxdx(1:3,3) = dSPHERICALdz
end subroutine get_jacobian
end module metric
