module io
   implicit none
   integer, parameter :: id = 1, nprocs = 1
end module io

module part
   implicit none
contains
   logical function isdead(i)
      integer, intent(in) :: i
      integer :: temp
      isdead = .false.
      temp = i
   end function isdead
end module part

module mpiutils
   implicit none
contains
   subroutine waitmyturn(id)
      integer, intent(in) :: id
      integer :: temp
      temp = id
   end subroutine waitmyturn
   
   subroutine endmyturn(id)
      integer, intent(in) :: id
      integer :: temp
      temp = id
   end subroutine endmyturn
end module mpiutils