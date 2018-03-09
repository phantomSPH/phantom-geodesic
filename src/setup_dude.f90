module setup

 implicit none

contains

!--- Setup a stick-figure person
subroutine setpart(xall,vall,np)
 use utils, only: cross_product, get_rotation_matrix
 integer, intent(out) :: np
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 integer :: i,index,nleg,nbody,narm,nhead
 real :: rotate_z(3,3),head_radius,rtan(3),r
 real, parameter :: pi = acos(-1.)
 real, parameter :: spacing=0.05, theta_z=pi/2.
 real, parameter :: translate(3) = (/20.,0.,0./)

 nbody = 10
 narm  = 5
 nhead = 10
 nleg=5
 np = nbody + (2*narm + 1) + nhead + 2*nleg
 head_radius=spacing*narm

 allocate(xall(3,np),vall(3,np))
 vall = 0.
 xall = 0.

 index = 0

 do i=1,nbody
    index = index+1
    xall(1,index)=spacing*(i-1)
 enddo

 do i=-narm,narm
    index = index+1
    xall(:,index) = (/0.75*(spacing*nbody),spacing*i,0./)
 enddo

 do i=1,nhead
    index = index+1
    xall(:,index) = (/nbody*spacing + head_radius + head_radius*cos(i*2.*pi/10.),head_radius*sin(i*2.*pi/10.),0./)
 enddo

 do i=1,nleg
    index=index+1
    xall(1,index) = -nleg*spacing+(i-1)*spacing
    xall(2,index) = -xall(1,index)
 enddo
 do i=1,nleg
    index=index+1
    xall(1,index) = -nleg*spacing+(i-1)*spacing
    xall(2,index) = xall(1,index)
 enddo

 call get_rotation_matrix(theta_z,rotate_z,'z')

 do i=1,np
    xall(:,i)=matmul(rotate_z,xall(:,i))
    xall(:,i)=xall(:,i) + translate
    call cross_product((/0.,0.,1./),xall(:,i),rtan)
    rtan = rtan/sqrt(dot_product(rtan,rtan)) ! Unit vector tangential to motion
    r = sqrt(dot_product(xall(:,i),xall(:,i)))
    ! vall(:,i) = rtan/sqrt(r)
 enddo

 ! vall(2,:) = 0.0521157

end subroutine setpart

end module setup
