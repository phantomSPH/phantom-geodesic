!----------------------------------------------------------------
!+
!  Evolve the position of a point particle in a metric (e.g. Schwarzschild)
!+
!----------------------------------------------------------------
program test
   use init, only: setup
   use metric, only: get_metric, mass1
   use force_gr, only: get_sourceterms
   use step, only: step_leapfrog, step_1, step_heuns
   use utils_gr, only: get_u0
   use output, only: write_out, write_ev, write_xyz, write_vxyz
   use checks, only: check
   implicit none
   integer :: np
   real, allocatable, dimension(:,:) :: xall,vall
   real, dimension(3):: x,v,fterm
   real :: time, energy_init, angmom_init, energy, angmom, U0, r
   real, parameter :: dt = 1.e-3, tmax = 10000., dtout = 100., dtout_ev = 10.
   integer :: nsteps, i,j, dnout, dnout_ev
   logical :: passed

   nsteps = int(tmax/dt)
   print*,'dt     = ',dt
   print*,'nsteps = ',nsteps
   time  = 0.
   dnout = int(dtout/dt)
   dnout_ev = int(dtout_ev/dt)
   print*,'START'
   call setup(xall, vall,np)

   angmom = 0.
   energy = 0.
   do i=1,np
      x = xall(:,i)
      v = vall(:,i)
      call check(x,v,passed)
      if (.not. passed) then
         STOP "Bad initial conditions!"
      endif

      ! For Schwarzschild only
      call get_u0(x,v,U0)
      r           = sqrt(dot_product(x,x))
      energy_init = energy_init + (1. - 2*mass1/r)*U0
      angmom_init = angmom_init + (x(1)*v(2)-x(2)*v(1))*U0
   enddo
   call write_out(time,xall,vall,np)

   do i=1,nsteps
      time = time + dt
      angmom =0.
      energy =0.

      do j=1,np

         x = xall(:,j)
         v = vall(:,j)
         call get_sourceterms(x,v,fterm)
         !call step_leapfrog(x,v,fterm,dt)
         call step_heuns(x,v,fterm,dt)
         !call step_1(x,v,fterm,dt)
         if (mod(i,dnout)==0) call check(x,v,passed)
         xall(:,j) = x
         vall(:,j) = v

         ! For Schwarzschild only
         call get_u0(x,v,U0)
         r           = sqrt(dot_product(x,x))
         energy = energy + (1. - 2*mass1/r)*U0
         angmom = angmom + (x(1)*v(2)-x(2)*v(1))*U0
      enddo

      if (mod(i,dnout_ev)==0) then
         call write_ev(time,energy,angmom)
         call write_xyz(time,xall,np)
         call write_vxyz(time,vall,np)
      endif
      if (mod(i,dnout)==0) then
         call check(x,v,passed)
         print*,'TIME =', time
         call write_out(time,xall,vall,np)
      endif

   enddo
end program test
