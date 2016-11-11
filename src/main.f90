!----------------------------------------------------------------
!+
!  Evolve the position of a point particle in a metric (e.g. Schwarzschild)
!+
!----------------------------------------------------------------
program test
   use init, only: setup,setup_dude,setup_sphere
   use metric, only: metric_type
   use force_gr, only: get_sourceterms
   use step, only: step_leapfrog, step_1, step_heuns, step_rk2
   use utils_gr, only: get_ev
   use output, only: write_out, write_ev, write_xyz, write_vxyz
   use checks, only: check,sanity_checks
   use utils_gr, only: get_rderivs
   implicit none
   integer :: np
   real, allocatable, dimension(:,:) :: xall,vall
   real, dimension(3):: x,v,fterm
   real :: time, energy_init, angmom_init, energy, angmom
   integer, parameter :: ndumps=1000!15000/15*20
   real, parameter :: dt = 1.e-4, tmax = 15000*2*20*0. + 310, dtout_ev = tmax/ndumps, dtout = dtout_ev*1000.
   integer :: nsteps, i,j, dnout, dnout_ev
   logical :: passed
   real :: start, finish, frac_done, twall_elapsed
   !real :: xblah(3),vblah(3),term1,term2,term3,term

   print*,"Metric type = ",metric_type
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
      call get_ev(x,v,energy,angmom)
      ! ! For Schwarzschild only
      ! call get_u0(x,v,U0)
      ! r           = sqrt(dot_product(x,x))
      ! energy = energy + (1. - rs/r)*U0
      ! angmom = angmom + (x(1)*v(2)-x(2)*v(1))*U0
      ! ! energy = energy + (1. - rs/x(1))*U0
      ! ! angmom = angmom + x(1)**2*v(3)*U0
   enddo
   energy_init = energy
   angmom_init = angmom

   call write_out(time,xall,vall,np)
   call write_ev(time,energy-energy_init,angmom-angmom_init)
   ! call write_ev(time,energy,angmom)
   call write_xyz(time,xall,np)
   call write_vxyz(time,vall,np)
   call cpu_time(start)
   do i=1,nsteps
      time = time + dt
      angmom =0.
      energy =0.
      do j=1,np
         x = xall(:,j)
         v = vall(:,j)
         call get_sourceterms(x,v,fterm)
         call step_leapfrog(x,v,fterm,dt)
         ! call step_heuns(x,v,fterm,dt)
         ! call step_rk2(x,v,fterm,dt)
         ! call step_1(x,v,fterm,dt)
         if (mod(i,dnout)==0) call check(x,v,passed)
         xall(:,j) = x
         vall(:,j) = v
         call get_ev(x,v,energy,angmom)
         ! ! For Schwarzschild only
         ! call get_u0(x,v,U0)
         ! r           = sqrt(dot_product(x,x))
         ! energy = energy + (1. - rs/r)*U0
         ! angmom = angmom + (x(1)*v(2)-x(2)*v(1))*U0
         ! ! energy = energy + (1. - rs/x(1))*U0
         ! ! angmom = angmom + x(1)**2*v(3)*U0
      enddo

      ! xblah = xall(:,1)
      ! vblah = vall(:,1)
      ! term1 = -rs/xblah(1)**2
      ! term2 = 2.*a*rs*vblah(3)/xblah(1)**2
      ! term3 = (2.*xblah(1)**5-a**2*xblah(1)**2*rs)*vblah(3)**2/xblah(1)**4
      ! term = (term1+term2+term3)*0.5*u0
      ! print*,'      fterm:',fterm(1),term,abs((term-fterm(1))/fterm(1))
      ! read*
      ! write(1,*) time,fterm

      if (mod(i,dnout_ev)==0) then
         call write_ev(time,energy-energy_init,angmom-angmom_init)
         ! call write_ev(time,energy,angmom)
         call write_xyz(time,xall,np)
         call write_vxyz(time,vall,np)
         twall_elapsed = finish-start
         frac_done = time/tmax
         print*,'t =', time,'%:',frac_done*100.,'t-minus (s):',(1.-frac_done)/frac_done*twall_elapsed
         call cpu_time(finish)
      endif
      if (mod(i,dnout)==0) then
         call check(x,v,passed)
         call write_out(time,xall,vall,np)
      endif
   enddo
end program test
