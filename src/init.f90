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

 call setpart(xall, vall,np)

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

end module init
