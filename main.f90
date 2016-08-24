!----------------------------------------------------------------
!+
!  Evolve the position of a point particle in a metric (e.g. Schwarzschild)
!+
!----------------------------------------------------------------
program test
 use init, only: initialise
 use metric, only: get_metric, mass1
 use force_gr, only: get_sourceterms
 use step, only: step_leapfrog, step_1, step_heuns
 use utils_gr, only: dot_product_gr
 use output, only: write_out
 use checks, only: check
 implicit none
 real, dimension(3):: x,v,fterm
 real, dimension(0:3,0:3) :: gcov, gcon
 real :: sqrtg, v4(0:3), time, energy_init, angmom_init, U0, r
 real, parameter :: dt = 1.e-3, tmax = 5000., dtout = 1.
 integer :: nsteps, i, dnout
 logical :: passed

 nsteps = int(tmax/dt)
 print*,'dt     = ',dt
 print*,'nsteps = ',nsteps
 time  = 0.
 dnout = int(dtout/dt)

 print*,'START'
 call initialise(x, v)
 call check(x,v,passed)
 if (.not. passed) then
  STOP "Bad initial conditions!"
 endif
 call get_metric(x,gcov,gcon,sqrtg)

 v4(0) = 1.
 v4(1:3) = v(1:3)
 U0 = 1./sqrt(-dot_product_gr(v4,v4,gcov))
 ! For Schwarzschild only
 r           = sqrt(dot_product(x,x))
 energy_init = (1. - 2*mass1/r)*U0
 angmom_init = (x(1)*v(2)-x(2)*v(1))*U0

 call write_out(time,x,v,energy_init,angmom_init)

 do i =1,nsteps
  time = time + dt
  call get_sourceterms(x,v,fterm)
  !call step_leapfrog(x,v,fterm,dt)
  call step_heuns(x,v,fterm,dt)
  !call step_1(x,v,fterm,dt)
  if (mod(i,dnout)==0) then
   call check(x,v,passed)
   print*,i, time
   call write_out(time,x,v,energy_init,angmom_init)
  endif
 enddo
end program test
