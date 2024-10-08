module set_geodesic
 implicit none

 real, parameter    :: pi = acos(-1.)
 integer, parameter :: ngtypes = 11
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
                       'parabola             ',&
                       'binary               '&
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
                       iparabola = 10,     &
                       ibinary  = 11
 real :: rp_newton, &
         inc_parabola
 integer :: gtype
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
subroutine setgeodesic(x,v,mall,np,type,r0)
 use metric,       only:metric_type, a, mass1 !,rs
 use metric_tools, only:coordinate_sys
 use force_gr,     only:get_sourceterms
 use utils_gr,     only:get_u0
 use utils,        only:get_rotation_matrix
 use prompting,    only:prompt
 real, intent(in), optional  :: r0
 integer, intent(in) :: type,np
 real, intent(in)    :: mall(np)
 real, intent(inout) :: x(3,np), v(3,np)
 real :: r, vy, x1
 real :: ra,va,omega,fac
 real :: rotate_y(3,3), inclination
 real :: theta,phi,m,q,rho2,y1,z1,vx,vz,rdot,thetadot
 real :: ecc,semia,rp,beta,rt
 real :: vhat(3),vmag,dx(3),dv(3),mtot
 real :: pos_mag,vel_mag,eng,semi_major,vel_orbit
 character(len=120)      :: filename
 integer                 :: ierr
 logical                 :: iexist

 print*,""
 print*,"We are using option: ",type

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
       x(1:3,np) = (/x1,0.,0./)
       v(1:3,np) = (/0.,vy,0./)
    case('Spherical')
       x(1:3,np) = (/r,0.5*pi,0./)
       v(1:3,np) = (/0.,0.,omega/)
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
       x(1:3,np) = (/x1,0.,0./)
    case('Spherical')
       x(1:3,np) = (/r,0.5*pi,0./)
    end select
    v(1:3,np) = (/0.,0.,0./)

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
       x(1:3,np) = (/ra,0.,0./)
       v(1:3,np) = (/0.,va,0./)
    case('Spherical')
       x(1:3,np) = (/ra,0.5*pi,0./)
       v(1:3,np) = (/0.,0.,va/ra /)
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
    x(1:3,np) = (/ra,0.,0./)
    v(1:3,np) = (/0.,va,0./)
    call get_rotation_matrix(inclination,rotate_y,'y')
    x(1:3,1) = matmul(rotate_y,x(1:3,1))

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
       x(1:3,np) = (/x1,0.,0./)
       v(1:3,np) = (/0.,fac*vy,0./)
    case('Spherical')
       x(1:3,np) = (/r,0.5*pi,0./)
       v(1:3,np) = (/0.,0.,fac*omega/)
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
       x(1:3,np) = (/x1,0.,0.+(fac-1.)/)
       v(1:3,np) = (/0.,vy,0./)
    case('Spherical')
       x(1:3,np) = (/r,0.5*pi*fac,0./)
       v(1:3,np) = (/0.,0.,omega/)
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
       x(1:3,np) = (/x1,y1,z1/)
       v(1:3,np) = (/vx,vy,vz/)
    case('Spherical')
       x(1:3,np) = (/r,theta,phi/)
       v(1:3,np) = (/rdot,thetadot,omega/)
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
      x(1:3,np) = (/x1,y1,z1/)
      v(1:3,np) = (/vx,vy,vz/)
    case('Spherical')
      STOP 'Need to be in cartesian'
    end select

 case(iellipse)
    ecc = 0.8
    rp  = 47.131 ! Tidal radius for solar type star around 1e6 Msun black hole
    inclination = 45.
    call prompt('eccentricity',ecc)
    !call prompt('r pericentre',rp)
    call prompt('semi-major', semia)
    call prompt('inclination (deg)',inclination)
    inclination = inclination/180. * pi
    !semia = rp/(1.-ecc)
    rp = semia*(1.-ecc)
    r  = semia*(1.+ecc)
    vy = sqrt(mass1*(1.-ecc)/r)
    print*,mass1,"mass1"
    x(1:3,np)  = (/r,0.,0./)
    v(1:3,np)  = (/0.,vy,0./)
    call get_rotation_matrix(-inclination,rotate_y,'y')
    x(1:3,1) = matmul(rotate_y,x(1:3,1))

    print*,'Period of orbit = ',2.*pi*sqrt(semia**3/1.)
    print*,'Suggested dt: ',(2.*pi*sqrt(rp**3))/100.
    print*,'Press ENTER to continue'
    read*

 case(iparabola)
    !default values
    rp_newton  = 47.131
    inc_parabola = 45.
    r = 500
    filename = 'orbit'//'.params'
    inquire(file=filename,exist=iexist)
    if (iexist) call read_setupfile(filename,ierr)
    if (.not. iexist .or. ierr /= 0) then
       call write_setupfile(filename)
       print*,' Edit '//trim(filename)//' and rerun'
       stop
    endif

    print*,rp_newton,"rp newton"
    y1 = -2.*rp_newton + r
    x1 = sqrt(r**2 - y1**2)
    x(1:3,np)  = (/x1,y1,0./)
    vmag = sqrt(2.*1./r)
    vhat = (/-2.*rp_newton,-x1,0./)/sqrt(4.*rp_newton**2 + x1**2)
    v(1:3,np)    = vmag*vhat
    inc_parabola = inc_parabola/180. * pi
    call get_rotation_matrix(-inc_parabola,rotate_y,'y')
    x(1:3,np) = matmul(rotate_y,x(1:3,np))
    v(1:3,np) = matmul(rotate_y,v(1:3,np))

   print*,'Suggested dt: ',(2.*pi*sqrt(rp_newton**3))/100.

 case(ibinary)
    mtot = sum(mall)
    call prompt('eccentricity',ecc)
    call prompt('semi-major',semia)

    ! set binary at apastron
    dx = (/semia*(1. + ecc),0.,0./)
    dv = (/0.,sqrt(semia*(1.-ecc**2)*mtot)/dx(1),0./)

    x(1:3,1) = -dx*mall(2)/mtot + (/10000000.,0.,0./)
    x(1:3,2) =  dx*mall(1)/mtot + (/10000000.,0.,0./)

    ! velocities
    v(1:3,1) = -dv*mall(2)/mtot
    v(1:3,2) =  dv*mall(1)/mtot

    ! The following parameters are a test collision case from Monte Carlo code of Alexander Heger
    ! x(1:3,1) = (/22170.50316413,  6506.44584756,    -0./)
    ! x(1:3,2) = (/22163.56404343,  6494.1490515 ,     0./)
    !
    ! v(1:3,1) = (/-0.00937932, -0.00118909, -0./)
    ! v(1:3,2) = (/-0.00903958, -0.00145593,  0./)

    ! this is the one we will test for different time step
    ! x(1:3,1) = (/ 4205172.308006175, 407876.95095306425,0./)
    ! x(1:3,2) = (/ 4205171.860429131, 407874.317039877, 0./)
    ! v(1:3,1) = (/  -0.02215735126131465, -0.0009752156643123136,0./)
    ! v(1:3,2) = (/ -0.021307350537264028, -0.001127753206150499,0./)

    ! x(1:3,1) = (/ 414635., 15098.168012089971,0./)
    ! x(1:3,2) = (/ 414634., 15100.455024693552, 0./)
    ! v(1:3,1) = (/ -0.0018125042852483984, 0.00017179484524338824,0./)
    ! v(1:3,2) = (/ -0.002577916936322562, -0.0002517095306686536,0./)

    ! x(1:3,1) = (/  6480927.435113909, 529582.2122989591,0./)
    ! x(1:3,2) = (/ 6480925.871844204, 529583.4294852862, 0./)
    ! v(1:3,1) = (/  -0.017214870274126038, -0.000318381598394309,0./)
    ! v(1:3,2) = (/ -0.017832134174461018, -0.0011111525441099524,0./)

    ! this one is inclined orbit.
    x(1:3,1) = (/4205171.490408985, 407876.40192287887, -0.9176798294906736/)
    v(1:3,1) = (/  -0.02188157081647087, -0.0007900244838645161, 0.0003095389204715863/)
    x(1:3,2) = (/ 4205172.678026321, 407874.8660700624, 0.9176798294906736/)
    v(1:3,2) = (/ -0.021583130982107807, -0.0013129443865982967, -0.0003095389204715863/)

 end select

 print*,""


end subroutine setgeodesic

!
!---Read/write setup file--------------------------------------------------
!
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20
 integer :: i

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# tde setup file'

 call write_inopt(rp_newton,'rp_newton','newtonian rp',iunit)
 call write_inopt(inc_parabola,'inc_parabola','inc of orbit',iunit)
 close(iunit)

end subroutine write_setupfile

subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",'reading setup options from '//trim(filename)
 nerr = 0
 ierr = 0
 call open_db_from_file(db,filename,iunit,ierr)

 call read_inopt(rp_newton,'rp_newton',db,min=0.,errcount=nerr)
 call read_inopt(inc_parabola,'inc_parabola',db,min=0.,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
     print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
     ierr = nerr
 endif
  print*,rp_newton,"rp_newton in readsetup"
end subroutine read_setupfile

end module set_geodesic
