program test
 use init, only: initialise
 use metric, only: get_sourceterms
 implicit none
 integer, parameter :: ndim = 3
 real(8), dimension(ndim):: x,v,pmom,fterm
 real(8), dimension(1+ndim,1+ndim) :: gcov, gcon
 !real(8) :: dt, tmax, time
 !integer :: nsteps, i
 !real(8) :: energy
 !real(8) :: angmom(ndim)

 print*,'hello world'

 call initialise(ndim, x, v)
 call get_sourceterms(ndim,x,v,gcov,gcon,fterm)

end program test
