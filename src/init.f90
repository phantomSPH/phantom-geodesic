module init
 use setup,    only:setpart
 use checks,   only:check
 use energies, only:get_ev
 implicit none

contains

subroutine initialise(xall,vall,np,energy,angmom)
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 integer, intent(out) :: np
 real,    intent(out) :: energy,angmom
 integer :: i
 real    :: x(3),v(3),energy_i,angmom_i
 logical :: passed

 if (command_argument_count()==0) then
    call setpart(xall,vall,np)
 else
    call read_phantom_ascii_dump(xall,vall,np)
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

end subroutine initialise

subroutine read_phantom_ascii_dump(xall,vall,np)
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 integer, intent(out) :: np
 character(len=32) :: start_dump
 character(len=1) :: hash
 integer :: i,nn(5)
 integer, parameter :: iunit=256
 real :: pmass,hsmooth,dens

 call get_command_argument(1,start_dump)
 print*,'Trying to start from dump file: ',start_dump
 open(unit=iunit,file=start_dump,status='old',action='read')

 nn = -1
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
! Particle positions and velocities
!
 do i = 1,np
    read(iunit,*) xall(1:3,i),pmass,hsmooth,dens,vall(1:3,i)
 enddo

 close(iunit)


end subroutine read_phantom_ascii_dump

end module init
