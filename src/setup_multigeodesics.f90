module setup
 use set_geodesic, only:setgeodesic

 implicit none

contains

!--- Setup up multiple test particles

subroutine setpart(xall,vall,np)
 integer, intent(out) :: np
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 real :: r0
 integer :: i
 np = 10
 allocate(xall(3,np),vall(3,np))
 do i=1,np
    r0 = 0.+i*4.
    call setgeodesic(xall(:,i),vall(:,i),'radial',r0)
 enddo
end subroutine setpart

end module setup
