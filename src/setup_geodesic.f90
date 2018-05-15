module setup
 use set_geodesic, only:setgeodesic,gtypelist,ngtypes,print_geodesic_choices,iprec
 use prompting,    only:prompt

 implicit none

contains

!--- Setup up a single test particle

subroutine setpart(xall,vall,np)
 integer, intent(out) :: np
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 integer :: gtype
 np = 1
 allocate(xall(3,np),vall(3,np))

 call print_geodesic_choices
 gtype = iprec
 call prompt(' Enter geodesic choice:',gtype)
 call setgeodesic(xall(:,1),vall(:,1),gtype)
end subroutine setpart

end module setup
