program test
 use init, only: initialise
 use metric, only: get_sourceterms
 implicit none
 !integer, parameter :: ndim = 3
 real, dimension(3):: x,v,pmom,fterm
 real, dimension(0:3,0:3) :: gcov, gcon
 !real(8) :: dt, tmax, time
 !integer :: nsteps, i
 !real(8) :: energy
 !real(8) :: angmom(ndim)

 print*,'hello world'

 call initialise(x, v)
 call get_sourceterms(x,v,fterm)

end program test
