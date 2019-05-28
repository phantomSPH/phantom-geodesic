module init
 use setup,    only:setpart
 use checks,   only:check
 use energies, only:get_ev
 implicit none

contains

subroutine initialise(time,xall,vall,np,energy,angmom)
 use output, only:write_out
 real, intent(in) :: time
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 integer, intent(out) :: np
 real,    intent(out) :: energy,angmom
 integer :: i
 real    :: x(3),v(3),energy_i,angmom_i
 logical :: passed

 if (command_argument_count()==0) then
    call setpart(xall,vall,np)
 else
    call read_dump(xall,vall,np)
 endif

 ! Check setup and compute total energy and angular momentum

 angmom = 0.
 energy = 0.
 do i=1,np
    x = xall(:,i)
    v = vall(:,i)
    call check(x,v,passed)
    if (.not. passed) then
       STOP "Bad initial conditions!"
    endif
    call get_ev(x,v,energy_i,angmom_i)
    energy = energy + energy_i
    angmom = angmom + angmom_i
 enddo
 print*,'Setup OK'
 print*,''

 if (command_argument_count()==0) then
    print*,'Writing to starting dump: output_00000.dat'
    call write_out(time,xall,vall,np)
    print*,'edit grtest.in file and run ./grtest output_00000.dat'
    STOP
 endif

end subroutine initialise

subroutine read_grtest_dump(iunit,xall,vall,np)
 use metric_tools, only: coordinate_sys
 use metric,       only: cartesian2spherical
 integer, intent(in) :: iunit
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 integer, intent(out) :: np
 character(len=1) :: hash
 integer :: i
 real :: x(3),v(3)

 np = -1

!
! Header tag and time
!
 do i=1,2
    read(iunit,*)
 enddo

!
! Number of particles
!
 read(iunit,'(a,i20)') hash,np
 allocate(xall(3,np),vall(3,np))

!
! Column labels
!
 read (iunit,*)

!
! Particle positions and velocities (assumed to be in cartesian coordinates)
!
 do i = 1,np
    read(iunit,*) x,v
    if (coordinate_sys == 'Spherical') then
       call cartesian2spherical(x,xall(:,i),v,vall(:,i))
    else
       xall(1:3,i) = x
       vall(1:3,i) = v
    endif
 enddo

end subroutine read_grtest_dump

subroutine read_phantom_ascii_dump(iunit,xall,vall,np)
 use metric_tools, only: coordinate_sys
 use metric,       only: cartesian2spherical
 integer, intent(in) :: iunit
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 integer, intent(out) :: np
 character(len=1) :: hash
 integer :: i
 real :: x(3),pmass,hsmooth,dens,v(3)

 np = -1

!
! Junk
!
 do i = 1,7
    read(iunit,*)
 enddo

!
! Number of particles
!
 read(iunit,'(a,i20)') hash,np
 allocate(xall(3,np),vall(3,np))

!
! Junk
!
 do i = 9,14
    read (iunit,*)
 enddo

!
! Particle positions and velocities (assumed to be in cartesian coordinates)
!
 do i = 1,np
    read(iunit,*) x,pmass,hsmooth,dens,v
    if (coordinate_sys == 'Spherical') then
      call cartesian2spherical(x,xall(:,i),v,vall(:,i))
    else
      xall(1:3,i) = x
      vall(1:3,i) = v
    endif
 enddo

end subroutine read_phantom_ascii_dump

subroutine read_dump(xall,vall,np)
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 integer, intent(out) :: np
 character(len=32)    :: tag
 character(len=32)    :: start_dump
 integer, parameter   :: iunit=256

 call get_command_argument(1,start_dump)
 print*,'Trying to start from dump file: ',start_dump

 open(unit=iunit,file=start_dump,status='old',action='read')
 read(iunit,'(a)') tag
 close(iunit)

 open(unit=iunit,file=start_dump,status='old',action='read')

 if (tag=='# grtest dumpfile') then
    call read_grtest_dump(iunit,xall,vall,np)
 else
    call read_phantom_ascii_dump(iunit,xall,vall,np)
 endif

 close(iunit)

end subroutine read_dump

end module init
