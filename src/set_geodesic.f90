module set_geodesic
 implicit none

 real, parameter :: pi = acos(-1.)

contains

!--- Several types of single particle geodesics
subroutine setgeodesic(x,v,type,r0)
 use metric, only: metric_type, a!,rs
 use metric_tools, only: coordinate_sys
 use force_gr, only: get_sourceterms
 use utils_gr, only: get_u0
 use utils,    only: get_rotation_matrix
 real, intent(in), optional  :: r0
 real, intent(out) :: x(3), v(3)
 real :: r, vy, x1
 real :: ra,va,omega,fac
 character(len=*), intent(in) :: type
 real :: rotate_y(3,3), inclination

 print*,""

 if (present(r0)) then
    r = r0
    write(*,'(a,f6.2)') ' Using init with r = ',r0
 endif

 select case(type)

 case('circular')
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

 case ('radial') ! Radial infall
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

 case('precession') ! Clement's orbit
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

 case('precession inclined')
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

 case('epicycle')
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

 case('vertical-oscillation')
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

 end select

 print*,""

end subroutine setgeodesic

end module set_geodesic