module setup
 use set_sphere, only:setsphere
 implicit none

contains

!--- Setup a proxy TDE via a sphere of particles

subroutine setpart(xall,vall,np)
 use utils, only: get_rotation_matrix
 integer, intent(out) :: np
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 real, parameter :: translate(3) = (/90.,0.,0./)
 real, parameter :: pi = acos(-1.)
 real    :: rotate_y(3,3)
 integer :: i

 call setsphere(xall,vall,np)

 call get_rotation_matrix(-pi/4.,rotate_y,'y')
 do i=1,np
    xall(:,i)=xall(:,i) + translate
    xall(:,i)=matmul(rotate_y,xall(:,i))
    vall(:,i)=(/0.,0.0521157,0./)
 enddo
end subroutine setpart

end module setup
