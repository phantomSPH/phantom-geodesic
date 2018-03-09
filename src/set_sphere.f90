module set_sphere

 implicit none

contains

!--- Setup a "sphere" of particles

subroutine setsphere(xall,vall,np)
 use utils, only: get_rotation_matrix
 integer, intent(out) :: np
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 real :: dlayer, dr
 integer :: i,j,k, nlayers, nrings,index,n,nringsmax
 real, parameter :: pi = acos(-1.)
 dlayer = 0.5
 dr = 0.5
 nlayers = 7
 nringsmax = 4
 index = 0
 np=107
 allocate(xall(3,np),vall(3,np))
 do i=0,(nlayers-1)/2 !layer index
    nrings = nringsmax-i
    do j=1,nrings !number of rings in layer
       n=5*(j-1)+1
       do k=1,n !number of points in ring
          index=index+1
          if (j==0) then
             xall(:,index)= (/0.,0.,i*dlayer/)
          else
             xall(:,index)= (/(j-1)*dr*cos(2.*pi/n*k),(j-1)*dr*sin(2.*pi/n*k), i*dlayer/)

             index=index+1
             xall(:,index)= (/(j-1)*dr*cos(2.*pi/n*k),(j-1)*dr*sin(2.*pi/n*k),-i*dlayer/)

          endif
       enddo
    enddo
 enddo
 vall = 0.
end subroutine setsphere

end module set_sphere
