module setup
 use set_geodesic, only:setgeodesic,iprec
 use prompting,    only:prompt

 implicit none

contains

!--- Setup up multiple test particles

subroutine setpart(xall,vall,np)
 integer, intent(out) :: np
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 real :: r0,r1,dr
 integer :: i, gtype

 !
 ! Defaults
 !
 np = 10
 r1 = 4.
 dr = 4.
 gtype = iprec

 call prompt(' Enter number of particles',np,0)
 call prompt(' Enter starting r1',r1,0.)
 call prompt(' Enter starting dr',dr,0.)
 call prompt(' Enter geodesic choice:',gtype)

 allocate(xall(3,np),vall(3,np))

 do i=1,np
    r0 = r1+(i-1)*dr
    call setgeodesic(xall(:,i),vall(:,i),gtype,r0)
 enddo

end subroutine setpart

end module setup
