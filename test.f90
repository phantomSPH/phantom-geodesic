program test
 use init, only: initialise
 use metric, only: get_sourceterms
 implicit none
 integer, parameter :: ndim = 3
 real, dimension(1:ndim):: x,v,pmom,fterm
 real, dimension(0:ndim,0:ndim) :: gcov, gcon
 !real(8) :: dt, tmax, time
 !integer :: nsteps, i
 !real(8) :: energy
 !real(8) :: angmom(ndim)

 print*,'hello world'

 call initialise(ndim, x, v)
 call get_sourceterms(ndim,x,v,fterm)

end program test
