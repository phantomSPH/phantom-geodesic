module setup
 use set_geodesic, only:setgeodesic

 implicit none

contains

!--- Setup up a single test particle

subroutine setpart(xall,vall,np)
 integer, intent(out) :: np
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 np = 1
 allocate(xall(3,np),vall(3,np))
 call setgeodesic(xall(:,1),vall(:,1),'precession')
end subroutine setpart

end module setup
