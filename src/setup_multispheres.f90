module setup
 use set_sphere, only:setsphere
 implicit none

contains

!--- Setup multiple spheres of particles

subroutine setpart(xall,vall,np)
 use utils,      only: get_rotation_matrix
 use set_sphere, only:setsphere
 integer, intent(out) :: np
 real, allocatable, intent(inout), dimension(:,:) :: xall,vall
 real, allocatable, dimension(:,:) :: x1,v1,x2,v2,x3,v3
 real, parameter :: pi = acos(-1.)
 real, parameter :: translate(3) = (/90.,0.,0./)
 integer :: n1, n2,i,n3
 real :: rotate_y(3,3), rotate_x(3,3), rotate_z(3,3),theta_y, theta_z, theta_x

 call setsphere(x1,v1,n1)
 call setsphere(x2,v2,n2)
 call setsphere(x3,v3,n3)
 np = n1+n2+n3

 call get_rotation_matrix(-pi/4.,rotate_y,'y')
 do i=1,n1
    x1(:,i)=x1(:,i) + translate
    x1(:,i)=matmul(rotate_y,x1(:,i))
    v1(:,i)=(/0.,0.0521157,0./)
 enddo
 do i=1,n2
    x2(:,i)=x2(:,i) + translate
    x2(:,i)=matmul(rotate_y,x2(:,i))
    v2(:,i)=(/0.,0.0521157,0./)
 enddo
 do i=1,n3
    x3(:,i)=x3(:,i) + translate
    x3(:,i)=matmul(rotate_y,x3(:,i))
    v3(:,i)=(/0.,0.0521157,0./)
 enddo

 theta_y=-0.5*pi
 call get_rotation_matrix(theta_y,rotate_y,'y')

 theta_z=-0.5*pi
 call get_rotation_matrix(theta_z,rotate_z,'z')

 theta_x= 0.5*pi
 call get_rotation_matrix(theta_x,rotate_x,'x')

 allocate(xall(3,np),vall(3,np))
 do i=1,n2
    x2(:,i)=matmul(rotate_y,x2(:,i))
    v2(:,i)=matmul(rotate_y,v2(:,i))
 enddo
 do i=1,n3
    x3(:,i)=matmul(rotate_x,matmul(rotate_z,x3(:,i)))
    v3(:,i)=matmul(rotate_x,matmul(rotate_z,v3(:,i)))
 enddo
 xall(:,1:n1) = x1
 vall(:,1:n1) = v1
 xall(:,n1+1:n1+n2) = x2
 vall(:,n1+1:n1+n2) = v2
 xall(:,n1+n2+1:n1+n2+n3) = x3
 vall(:,n1+n2+1:n1+n2+n3) = v3

end subroutine setpart

end module setup
