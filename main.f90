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
 real :: sqrtg, v4(0:3), time, energy_init, angmom_init, U0
 real, parameter :: dt = 1.e-3, tmax = 5000., dt_out=.1
 integer :: nsteps, i, n_out
 nsteps = int(tmax/dt)
 print*,'dt = ',dt
 print*,'nsteps = ',nsteps
 time = 0.
 n_out=int(dt_out/dt)
 !call sanity_checks
 !stop "sanity checks"

 print*,'START'
 call initialise(x, v)

 !Check 'causality condition thingamajig'
 call get_metric(x,gcov,gcon,sqrtg)
 v4(0) = 1.
 v4(1:3) = v(1:3)
 if (dot_product_gr(v4,v4,gcov)>0.) then
  print*,"dot_product_gr(v4,v4,gcov)=",dot_product_gr(v4,v4,gcov)
  STOP "Bad initial conditions!"
 endif
 U0 = 1./sqrt(-dot_product_gr(v4,v4,gcov))
 
 ! For Schwarzschild only
 energy_init = (1. - 2*mass1/sqrt(dot_product(x,x)))/sqrt(-dot_product_gr(v4,v4,gcov))
 angmom_init = (x(1)*v(2)-x(2)*v(1))*U0
 
 call write_out(time,x,v,energy_init,angmom_init)

 do i =1,nsteps
  time = time + dt
  call check(x,v)
  call get_sourceterms(x,v,fterm)
  !call step_leapfrog(x,v,fterm,dt)
  call step_heuns(x,v,fterm,dt)
  !call step_1(x,v,fterm,dt)
  if (mod(i,n_out)==0) then
   print*,i, time
   !call check(x,v)
   call write_out(time,x,v,energy_init,angmom_init)
  endif
 enddo
end program test
