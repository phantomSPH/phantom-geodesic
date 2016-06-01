program test
 use init, only: initialise
 use metric, only: get_sourceterms
 use step, only: step_leapfrog, step_1, step_heuns
 implicit none
 !integer, parameter :: ndim = 3
 real, dimension(3):: x,v,pmom,fterm
 real, dimension(0:3,0:3) :: gcov, gcon
 real(8) :: dt, tmax, time
 integer :: nsteps, i
 !real(8) :: energy
 !real(8) :: angmom(ndim)
 nsteps = 500

 !print*,'hello world'

 call initialise(x, v)
 call get_sourceterms(x,v,fterm)
 !print*, x, v, fterm

 dt = 1.e-7
 do i =1,nsteps
  print*,i
  !call step_leapfrog(3,x,v,fterm,dt)
  !call step_heuns(x,v,fterm,dt)
  call step_1(x,v,fterm,dt)
  write(1,*) x
  !print*, x,v,fterm
 enddo
end program test
