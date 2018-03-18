module setup
 use set_geodesic, only:setgeodesic
 use prompting,    only:prompt

 implicit none

contains

!--- Setup up a single test particle

subroutine setpart(xall,vall,np)
 integer, intent(out) :: np
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 character(len=20) :: gtype
 np = 1
 allocate(xall(3,np),vall(3,np))

 gtype = 'precession'
 call prompt(' Enter geodesic choice:',gtype)
 call setgeodesic(xall(:,1),vall(:,1),trim(gtype))
end subroutine setpart

end module setup
