module utils
 implicit none

contains

!--- Subroutine to compute the 3-vector cross product
subroutine cross_product(a,b,c)
 real, dimension(3), intent(in) :: a,b
 real, dimension(3), intent(out) :: c

 c(1) = a(2)*b(3) - a(3)*b(2)
 c(2) = a(3)*b(1) - a(1)*b(3)
 c(3) = a(1)*b(2) - a(2)*b(1)

end subroutine cross_product

!--- Subroutine to return a rotation matrix, given some angle, about
!    one of the cartesian axes.
subroutine get_rotation_matrix(angle,rotation_matrix,axis)
 real, intent(in)  :: angle
 character(len=*), intent(in) :: axis
 real, intent(out) ::  rotation_matrix(3,3)
 real, dimension(3,3) :: rotate_x,rotate_y,rotate_z

 select case(axis)
 case('x')
    rotate_x(1,:) = (/1. , 0.        , 0.        /)
    rotate_x(2,:) = (/0. , cos(angle),-sin(angle)/)
    rotate_x(3,:) = (/0. , sin(angle), cos(angle)/)
    rotation_matrix = rotate_x

 case('y')
    rotate_y(1,:) = (/ cos(angle), 0. , sin(angle)/)
    rotate_y(2,:) = (/ 0.        , 1. , 0.        /)
    rotate_y(3,:) = (/-sin(angle), 0. , cos(angle)/)
    rotation_matrix = rotate_y

 case('z')
    rotate_z(1,:) = (/cos(angle),-sin(angle),0./)
    rotate_z(2,:) = (/sin(angle), cos(angle),0./)
    rotate_z(3,:) = (/0.        ,0.         ,1./)
    rotation_matrix = rotate_z
 end select

end subroutine get_rotation_matrix

subroutine timer(time)
 real, intent(out) :: time
 integer :: cc,cr,cm

 call system_clock(count=cc,count_rate=cr,count_max=cm)
 time = real(cc)/real(cr)

end subroutine timer

end module utils
