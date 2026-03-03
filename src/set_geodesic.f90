module set_geodesic
 implicit none

 real, parameter    :: pi = acos(-1.)
 integer, parameter :: ngtypes = 13
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
                       'binary               ',&
                       'single               ',&
                       'emilio               '&
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
                       ibinary  = 11,      &
                       isingle  = 12,      &
                       iemilio  = 13

 real :: rp_newton, &
         inc_parabola, &
         beta, &
         phi
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
 real :: r, vy, x1, x_shift(3), v_shift(3)
 real :: ra,va,omega,fac,dv_mag
 real :: rotate_y(3,3), inclination
 real :: theta,m,q,rho2,y1,z1,vx,vz,rdot,thetadot
 real :: ecc,semia,rp,rt
 real :: vhat(3),vmag,dx(3),dv(3),mtot
 real :: xcm(3),vcm(3),phidot,angmom,mb,ro
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
    rp_newton  = 47.151
    inc_parabola = 45.

    filename = 'orbit'//'.params'
    inquire(file=filename,exist=iexist)
    if (iexist) call read_setupfile(filename,ierr)
    if (.not. iexist .or. ierr /= 0) then
       call write_setupfile(filename)
       print*,' Edit '//trim(filename)//' and rerun'
       stop
    endif
    r = 1e5 ! start far away, so that we can set a parabolic orbit with the Newtonian rp = rp_newton in GR

    y1 = -2.*rp_newton + r
    x1 = sqrt(r**2 - y1**2)
    x(1:3,np)  = (/x1,y1,0./)
    vmag = sqrt(2.*mass1/r)
    vhat = (/-2.*rp_newton,-x1,0./)/sqrt(4.*rp_newton**2 + x1**2)
    v(1:3,np)    = vmag*vhat
    inc_parabola = inc_parabola/180. * pi
    call get_rotation_matrix(-inc_parabola,rotate_y,'y')
    x(1:3,np) = matmul(rotate_y,x(1:3,np))
    v(1:3,np) = matmul(rotate_y,v(1:3,np))

    print*, x(1:3,np), "x(1:3,np)", v(1:3,np), "v(1:3,np)"

   print*,'Suggested dt: ',(2.*pi*sqrt(rp_newton**3))/100.

 case(ibinary)
    mtot = sum(mall)

    semia = 4.71307 ! test value
    ecc = 0.00

    call prompt('eccentricity',ecc)
    call prompt('semi-major',semia)

    beta = 1.0

    rt = semia * (mass1 / mtot)**(1./3.)
    rp = rt / beta
    r = 1000000.

    y1 = -2.*rp + r
    x1 = sqrt(r**2 - y1**2)

    x_shift(1:3)  = (/x1,y1,0./)
    vmag = sqrt(2.*1./r)

    vhat = (/-2.*rp,-x1,0./)/sqrt(4.*rp**2 + x1**2)
    v_shift(1:3)  = vmag*vhat

    ! r_apo = a * (1 + e)
    dx = (/ 0., semia * (1.0 + ecc), 0. /)

    dv_mag = sqrt( mtot * (1.0 - ecc) / (semia * (1.0 + ecc)) )
    dv = (/ dv_mag, 0., 0. /)

    if (mass1 == 0.) then
      x(1:3,1) = -dx * mall(2) / mtot
      x(1:3,2) =  dx * mall(1) / mtot

      v(1:3,1) = -dv * mall(2) / mtot + (/0.0, 0.5, 0.0/) ! add a small velocity to have a drift
      v(1:3,2) =  dv * mall(1) / mtot + (/0.0, 0.5, 0.0/)
    else 
      x(1:3,1) = -dx * mall(2) / mtot + x_shift
      x(1:3,2) =  dx * mall(1) / mtot + x_shift    

      v(1:3,1) = -dv * mall(2) / mtot + v_shift
      v(1:3,2) =  dv * mall(1) / mtot + v_shift
    endif

  
 case(isingle)
   ! can provide the initial position and velocity of a single particle, and it will just follow the geodesic. Useful for testing.
   print*,'#--- Single particle ---#'
   x(1:3,1) = (/1466162.2748927348, 319546.2615316531, -0.0/)
   v(1:3,1) = (/-0.035833985056138815, -0.0036875703926903193, -0.0/)


 case(iemilio)

    filename = 'sample'//'.params'
    inquire(file=filename,exist=iexist)
    if (iexist) call read_setupfile_sampling(filename,ierr)
    if (.not. iexist .or. ierr /= 0) then
       call write_setupfile_sampling(filename)
       print*,' Edit '//trim(filename)//' and rerun'
       stop
    endif

   a = 1.0130607675019032 ! 0.01 AU in code units
   mb = mall(1) + mall(2)
   rt = a*(mass1/mb)**(1./3.)
   phidot = - sqrt(mb/a**3) ! considering retrograde orbit
   rp = rt/beta
   angmom = sqrt(2*mass1*rp**2 / (rp - 2*mass1))
   ro = 50*rt

   print*, beta, 'beta', rp, 'rp', rt, 'rt'

   xcm = (/ro,0.,0./)
   ! using the centre of mass velocity in Schwarzschild metric
   vcm = (/-(1-(2*mass1/ro))*sqrt(2*mass1/ro - (angmom**2/ro**2)*(1-(2*mass1/ro))), (1-(2*mass1/ro))*angmom/ro, 0./)

   x(1,1) = xcm(1) + mall(1)/mb * a * cos(phi)
   x(2,1) = xcm(2) + mall(1)/mb * a * sin(phi)
   x(3,1) = 0.
   
   x(1,2) = xcm(1) - mall(2)/mb * a * cos(phi)
   x(2,2) = xcm(2) - mall(2)/mb * a * sin(phi)
   x(3,2) = 0.
   
   v(1,1) = vcm(1) - mall(1)/mb * a * phidot * sin(phi)
   v(2,1) = vcm(2) + mall(1)/mb * a * phidot * cos(phi)
   v(3,1) = 0.
   
   v(1,2) = vcm(1) + mall(2)/mb * a * phidot * sin(phi)
   v(2,2) = vcm(2) - mall(2)/mb * a * phidot * cos(phi)
   v(3,2) = 0.

   
   print*,x(1:3,1),'1st star'
   print*,x(1:3,2),'2nd star'
   print*,v(1:3,1),'1st vel'
   print*,v(1:3,2),'2nd vel'


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

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# tde setup file'

 call write_inopt(rp_newton,'rp_newton','newtonian rp',iunit)
 call write_inopt(inc_parabola,'inc_parabola','inc of orbit',iunit)
 close(iunit)

end subroutine write_setupfile

subroutine write_setupfile_sampling(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# tde setup file'

 call write_inopt(beta,'beta','beta',iunit)
 call write_inopt(phi,'phi','phi',iunit)
 close(iunit)

end subroutine write_setupfile_sampling

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

end subroutine read_setupfile

subroutine read_setupfile_sampling(filename,ierr)
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

 call read_inopt(beta,'beta',db,min=0.,errcount=nerr)
 call read_inopt(phi,'phi',db,min=0.,errcount=nerr)
 call close_db(db)
 if (nerr > 0) then
     print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
     ierr = nerr
 endif

end subroutine read_setupfile_sampling

end module set_geodesic
