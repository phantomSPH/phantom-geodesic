module init
 use setup,    only:setpart,get_binary_mass
 use checks,   only:check
 use energies, only:get_ev,get_newtonian_energy
 implicit none

contains

subroutine initialise(time,xall,vall,np,energy,angmom,mall)
 use output, only:write_out
 real, intent(in) :: time
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 real, allocatable, intent(inout), dimension(:)   :: mall
 integer, intent(out) :: np
 real,    intent(out) :: energy,angmom
 integer :: i
 real    :: x(3),v(3),energy_i,angmom_i
 logical :: passed


 if (command_argument_count()==0) then
    call setpart(xall,vall,np,mall)
    if (np == 1) then
      allocate(mall(np))
      mall(:) = 0.
    endif
 else
    call read_dump(xall,vall,np,mall)
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
 call get_newtonian_energy(np,xall,vall,energy,mall)
 print*,'Setup OK'
 print*,''

 if (command_argument_count()==0) then
    print*,'Writing to starting dump: output_00000.dat'
    call write_out(time,xall,vall,np,mall)
    print*,'edit grtest.in file and run ./grtest output_00000.dat'
    STOP
 endif

end subroutine initialise

subroutine read_grtest_dump(iunit,xall,vall,np,mall)
 use metric_tools, only: coordinate_sys
 use metric,       only: cartesian2spherical
 integer, intent(in) :: iunit
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 real, allocatable, intent(inout), dimension(:) :: mall
 integer, intent(out) :: np
 character(len=1) :: hash
 integer :: i
 real :: x(3),v(3),m

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
 allocate(xall(3,np),vall(3,np),mall(np))

!
! Column labels
!
 read (iunit,*)

!
! Particle positions and velocities (assumed to be in cartesian coordinates)
!
 do i = 1,np
    read(iunit,*) x,v,m
    if (coordinate_sys == 'Spherical') then
       call cartesian2spherical(x,xall(:,i),v,vall(:,i))
    else
       xall(1:3,i) = x
       vall(1:3,i) = v
       mall(i) = m
    endif
 enddo

end subroutine read_grtest_dump

subroutine read_phantom_ascii_dump(iunit,xall,vall,np,mall)
 use metric_tools, only: coordinate_sys
 use metric,       only: cartesian2spherical
 integer, intent(in) :: iunit
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 real, allocatable, intent(inout), dimension(:) :: mall
 integer, intent(out) :: np
 character(len=1) :: hash
 integer :: i
 real :: x(3),pmass,hsmooth,dens,v(3),m

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
 allocate(xall(3,np),vall(3,np),mall(np))

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
      mall(i) = pmass
    endif
 enddo

end subroutine read_phantom_ascii_dump

subroutine read_dump(xall,vall,np,mall)
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 real, allocatable, intent(inout), dimension(:) :: mall
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
    call read_grtest_dump(iunit,xall,vall,np,mall)
 else
    call read_phantom_ascii_dump(iunit,xall,vall,np,mall)
 endif

 close(iunit)

end subroutine read_dump

end module init
