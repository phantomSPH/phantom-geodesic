module set_geodesic
 implicit none

 real, parameter    :: pi = acos(-1.)
 integer, parameter :: ngtypes = 10
 character(len=*), parameter  :: &
  gtypelist(ngtypes) = (/&
                       'circular             ',&
                       'radial               ',&
                       'precession           ',&
                       'precession inclined  ',&
                       'epicycle             ',&
                       'vertical-oscillation ',&
                       'circular-inclined    ',&
                       'custom               ',&
                       'ellipse              ',&
                       'parabola             ' &
                       /)

 integer, parameter :: &
                       icirc    = 1,       &
                       irad     = 2,       &
                       iprec    = 3,       &
                       iprecinc = 4,       &
                       iepi     = 5,       &
                       ivert    = 6,       &
                       icircinc = 7,       &
                       icustom  = 8,       &
                       iellipse = 9,       &
                       iparabola = 10

contains

subroutine print_geodesic_choices()
 integer :: i
 print*,''
 print*,'---------------------'
 print*,'Geodesic choices:'
 do i=1,ngtypes
    write(*,'(i2,")  ",a)') i,gtypelist(i)
 enddo
end subroutine print_geodesic_choices

!--- Several types of single particle geodesics
subroutine setgeodesic(x,v,type,r0)
 use metric,       only:metric_type, a!,rs
 use metric_tools, only:coordinate_sys
 use force_gr,     only:get_sourceterms
 use utils_gr,     only:get_u0
 use utils,        only:get_rotation_matrix
 use prompting,    only:prompt
 real, intent(in), optional  :: r0
 real, intent(out) :: x(3), v(3)
 integer, intent(in) :: type
 real :: r, vy, x1
 real :: ra,va,omega,fac
 real :: rotate_y(3,3), inclination
 real :: theta,phi,m,q,rho2,y1,z1,vx,vz,rdot,thetadot
 real :: ecc,semia,rp
 real :: vhat(3),vmag

 print*,""

 if (present(r0)) then
    r = r0
    write(*,'(a,f6.2)') ' Using init with r = ',r0
 endif

 select case(type)

 case(icirc)
    print*,'#--- Circular velocity in x-y plane, anticlockwise ---#'
    if (.not. present(r0)) then
       print*,'Enter radius r for (r,theta,phi)=(r,pi/2,0):'
       read(*,*) r
    endif
    if (metric_type=='Schwarzschild') then
       x1 = r                                    ! x1 = r in Schwarzschild
       omega = sqrt(1./r**3)
       vy = sqrt(1./r)
    elseif(metric_type=='Kerr') then
       ! r = x1 !sqrt(x1**2-a**2)                ! x1 /= r in Kerr
       omega = 1./(r**(1.5)+a)
       x1 = sqrt(r**2 + a**2)
       vy = x1*omega
    elseif(metric_type=='Minkowski') then
       STOP 'Cannot make circular orbits in Minkowski metric.'
    endif
    select case(coordinate_sys)
    case('Cartesian')
       x = (/x1,0.,0./)
       v = (/0.,vy,0./)
    case('Spherical')
       x = (/r,0.5*pi,0./)
       v = (/0.,0.,omega/)
    end select
    print*,'period =',2*pi/omega
    print*,'Press ENTER to continue:'
    read*

 case (irad) ! Radial infall
    print*,'#--- Radial infall ---#'
    if (.not. present(r0)) then
       print*,'Enter radius r for (r,theta,phi)=(r,pi/2,0):'
       read(*,*) r
    endif
    select case(coordinate_sys)
    case('Cartesian')
       x1 = sqrt(r**2 + a**2)
       x = (/x1,0.,0./)
    case('Spherical')
       x = (/r,0.5*pi,0./)
    end select
    v = (/0.,0.,0./)

 case(iprec) ! Clement's orbit
    ra = 90.
    va = 0.0521157 ! velocity giving a pericenter rp = 10
    print*,'#--- Precessing orbit: r =',ra,' and vy = ',va
    if (present(r0)) then
       print*, "WARNING: You're trying to set r = ",r0," but r is already set for this type of orbit. Continue?"
       read*,
    endif
    select case(coordinate_sys)
    case('Cartesian')
       x = (/ra,0.,0./)
       v = (/0.,va,0./)
    case('Spherical')
       x = (/ra,0.5*pi,0./)
       v = (/0.,0.,va/ra /)
       if (.not. metric_type=='Schwarzschild') STOP 'Only have precession setup for spherical in Schwarzschild'
    end select

 case(iprecinc)
    ra = 90.
    va = 0.0521157 ! velocity giving a pericenter rp = 10
    inclination = -pi/6
    print*,'#--- Inclined precessing orbit: r =',ra,' and vy = ',va,' inclined to plane by ',inclination*180./pi, 'degrees.'
    if (present(r0)) then
       print*, "WARNING: You're trying to set r = ",r0," but r is already set for this type of orbit. Continue?"
       read*,
    endif
    if (.not. metric_type=='Kerr') then
       print*,"Warning, you are not using the Kerr metric for 'precession inclined'...result will be same as Schwarzschild"
       read*
    endif
    if (.not. coordinate_sys=='Cartesian') STOP "Haven't tested 'precession inclined' for Spherical coordinates"
    x = (/ra,0.,0./)
    v = (/0.,va,0./)
    call get_rotation_matrix(inclination,rotate_y,'y')
    x = matmul(rotate_y,x)

 case(iepi)
    print*,'#--- Radial epicyclic motion in x-y plane, anticlockwise ---#'
    if (.not. present(r0)) then
       print*,'Enter radius r for (r,theta,phi)=(r,pi/2,0):'
       read(*,*) r
    endif
    if (metric_type=='Schwarzschild') then
       x1 = r                                    ! x1 = r in Schwarzschild
       omega = sqrt(1./r**3)
       ! vy    = sqrt(1./r)
       vy    = x1*omega
    elseif(metric_type=='Kerr') then
       ! r = x1 !sqrt(x1**2-a**2)                ! x1 /= r in Kerr
       omega = 1./(r**(1.5)+a)
       x1 = sqrt(r**2 + a**2)
       vy = x1*omega
    elseif(metric_type=='Minkowski') then
       STOP 'Cannot make circular orbits in Minkowski metric.'
    endif
    fac = 1.00001
    select case(coordinate_sys)
    case('Cartesian')
       x = (/x1,0.,0./)
       v = (/0.,fac*vy,0./)
    case('Spherical')
       x = (/r,0.5*pi,0./)
       v = (/0.,0.,fac*omega/)
    end select
    print*,'period =',2*pi/omega

 case(ivert)
    print*,'#--- Small vertical-oscillation from circular orbit in x-y plane, anticlockwise ---#'
    if (.not. present(r0)) then
       print*,'Enter radius r for (r,theta,phi)=(r,pi/2,0):'
       read(*,*) r
    endif
    if (metric_type=='Schwarzschild') then
       x1 = r                                    ! x1 = r in Schwarzschild
       omega = sqrt(1./r**3)
       ! vy    = sqrt(1./r)
       vy    = x1*omega
    elseif(metric_type=='Kerr') then
       ! r = x1 !sqrt(x1**2-a**2)                ! x1 /= r in Kerr
       omega = 1./(r**(1.5)+a)
       x1 = sqrt(r**2 + a**2)
       vy = x1*omega
    elseif(metric_type=='Minkowski') then
       STOP 'Cannot make circular orbits in Minkowski metric.'
    endif
    fac = 1.00001
    select case(coordinate_sys)
    case('Cartesian')
       x = (/x1,0.,0.+(fac-1.)/)
       v = (/0.,vy,0./)
    case('Spherical')
       x = (/r,0.5*pi*fac,0./)
       v = (/0.,0.,omega/)
    end select

 case(icircinc)
    print*,'#--- Circle incline to z=0 plane ---#'
    if(metric_type=='Minkowski') STOP 'Cannot make circular orbits in Minkowski metric.'
    if (.not. present(r0)) then
       r = 50.
       call prompt('Enter radius r for (r,theta,phi)=(r,theta,0)',r,0.)
    endif
    theta    = 45.
    call prompt('Enter theta (inclination in degrees from z-axis)',theta,0.,180.)
    theta    = theta*pi/180. ! convert to radians
    phi      = 0.
    m        = 1.
    q        = sqrt(r**2 - a**2*cos(theta)**2)
    rho2     = r**2 + a**2*cos(theta)**2
    omega    = q*sqrt(m)/(sin(theta)*(rho2*sqrt(r)+a*q*sqrt(m)*sin(theta))) !shakura 1987
    rdot     = 0.
    thetadot = 0.

    ! Cartesian coordinates
    x1 = sqrt(r**2+a**2)*sin(theta)*cos(phi)
    y1 = sqrt(r**2+a**2)*sin(theta)*sin(phi)
    z1 = r*cos(theta)
    vx = r/sqrt(r**2+a**2)*sin(theta)*cos(phi)*rdot + sqrt(r**2+a**2)*(cos(theta)*cos(phi)*thetadot-sin(theta)*sin(phi)*omega)
    vy = r/sqrt(r**2+a**2)*sin(theta)*sin(phi)*rdot + sqrt(r**2+a**2)*(cos(theta)*sin(phi)*thetadot+sin(theta)*cos(phi)*omega)
    vz = cos(theta)*rdot-r*sin(theta)*thetadot

    select case(coordinate_sys)
    case('Cartesian')
       x = (/x1,y1,z1/)
       v = (/vx,vy,vz/)
    case('Spherical')
       x = (/r,theta,phi/)
       v = (/rdot,thetadot,omega/)
    end select
    print*,'period =',2*pi/abs(omega)
    print*,'Press ENTER to continue:'
    read*

 case (icustom) ! custom setup
    x1 = 0.
    y1 = 0.
    z1 = 0.
    vx = 0.
    vy = 0.
    vz = 0.
    call prompt('x',x1)
    call prompt('y',y1)
    call prompt('z',z1)
    call prompt('vx',vx)
    call prompt('vy',vy)
    call prompt('vz',vz)
    select case(coordinate_sys)
    case('Cartesian')
      x = (/x1,y1,z1/)
      v = (/vx,vy,vz/)
    case('Spherical')
      STOP 'Need to be in cartesian'
    end select

 case(iellipse)
    ecc = 0.8
    rp  = 47.131 ! Tidal radius for solar type star around 1e6 Msun black hole
    inclination = 45.
    call prompt('eccentricity',ecc)
    call prompt('r pericentre',rp)
    call prompt('inclination (deg)',inclination)
    inclination = inclination/180. * pi
    semia = rp/(1.-ecc)
    r  = semia*(1.+ecc)
    vy = sqrt((1.-ecc)/r)
    x  = (/r,0.,0./)
    v  = (/0.,vy,0./)
    call get_rotation_matrix(-inclination,rotate_y,'y')
    x = matmul(rotate_y,x)

    print*,'Period of orbit = ',2.*pi*sqrt(semia**3/1.)
    print*,'Suggested dt: ',(2.*pi*sqrt(rp**3))/100.
    print*,'Press ENTER to continue'
    read*

 case(iparabola)
    rp  = 47.131 ! Tidal radius for solar type star around 1e6 Msun black hole
    r = 500. !2.*rp
    call prompt('r pericentre',rp)
    if (r < 2.*rp) then
       r = 2.*rp
    endif

    y1 = -2.*rp + r
    x1 = sqrt(r**2 - y1**2)
    x  = (/x1,y1,0./)
    vmag = sqrt(2.*1./r)
    vhat = (/-2.*rp,-x1,0./)/sqrt(4.*rp**2 + x1**2)
    v    = vmag*vhat

    inclination = 45.
    call prompt('inclination (deg)',inclination)
    inclination = inclination/180. * pi
    call get_rotation_matrix(-inclination,rotate_y,'y')
    x = matmul(rotate_y,x)
    v = matmul(rotate_y,v)

   print*,'Suggested dt: ',(2.*pi*sqrt(rp**3))/100.
   print*,'Press ENTER to continue'
   read*

 end select

 print*,""

end subroutine setgeodesic

end module set_geodesic
