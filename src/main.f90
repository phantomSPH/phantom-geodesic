!----------------------------------------------------------------
!+
!  Evolve the position of a point particle in a metric (e.g. Schwarzschild)
!+
!----------------------------------------------------------------
program test
use init,         only: setup
use metric,       only: metric_type
use metric_tools, only: coordinate_sys
use force_gr,     only: get_sourceterms
use step,         only: timestep, step_type
use utils_gr,     only: get_ev
use output,       only: write_out, write_ev, write_xyz, write_vxyz
use checks,       only: check,sanity_checks
use utils_gr,     only: get_rderivs
implicit none

integer, parameter :: ndumps=15000
real,    parameter :: dt = 1.e-2, tmax = 30000., dtout_ev = tmax/ndumps, dtout = dtout_ev*1000.

real, allocatable, dimension(:,:) :: xall,vall
real, dimension(3) :: x,v
integer :: np
real    :: time, energy_init, angmom_init, energy, angmom, energy_i, angmom_i
integer :: nsteps, i,j, dnout, dnout_ev
logical :: passed
real    :: start, finish, frac_done, twall_elapsed
real    :: tminus
integer :: percentage,prev_percent

   print*,'-------------------------------------------------------------------'
   print*,'GR-TEST'
   print*,'-------------------------------------------------------------------'
   print*,               'Metric type       = ',metric_type
   print*,               'Coord. sys. type  = ',coordinate_sys
   print*,               'Timestepping used = ',step_type
   write(*,'(a,f10.2)') ' dt                = ',dt
   write(*,'(a,f10.2)') ' tmax              = ',tmax
   write(*,'(a,f10.2)') ' dtout_ev          = ',dtout_ev
   nsteps   = int(tmax/dt)
   dnout    = int(dtout/dt)
   dnout_ev = int(dtout_ev/dt)
   print*,'-------------------------------------------------------------------'
   print*,'START'
   print*,'-------------------------------------------------------------------'
   time     = 0.
   call setup(xall, vall,np)
   print*,'Ready...'
   read*

   angmom = 0.
   energy = 0.
   do i=1,np
      x = xall(:,i)
      v = vall(:,i)
      call check(x,v,passed)
      if (.not. passed) then
         STOP "Bad initial conditions!"
      endif
      call get_ev(x,v,energy_i,angmom_i)
      energy = energy + energy_i
      angmom = angmom + angmom_i
   enddo
   energy_init = energy
   angmom_init = angmom

   call write_out(time,xall,np)
   call write_ev(time,energy-energy_init,angmom-angmom_init)
   call write_xyz(time,xall,np)
   call write_vxyz(time,vall,np)
   call cpu_time(start)
   do i=1,nsteps

      angmom =0.
      energy =0.
      do j=1,np
         x = xall(:,j)
         v = vall(:,j)
         call timestep(time,dt,x,v)
         if (mod(i,dnout)==0) call check(x,v,passed)
         xall(:,j) = x
         vall(:,j) = v
         call get_ev(x,v,energy_i,angmom_i)
         energy = energy + energy_i
         angmom = angmom + angmom_i
      enddo

      prev_percent = 0
      if (mod(i,dnout_ev)==0) then
         call write_ev(time,energy-energy_init,angmom-angmom_init)
         call write_xyz(time,xall,np)
         call write_vxyz(time,vall,np)
         twall_elapsed = finish-start
         frac_done     = time/tmax
         percentage    = nint(frac_done*100.)
         tminus        = (1.-frac_done)/frac_done*twall_elapsed
         if (percentage == prev_percent + 1) then
            write(*,'(i4.1,a,f10.2,a,f10.2)') percentage,'%   t =', time,'      t-minus (s):',tminus
         endif
         prev_percent = percentage
         call cpu_time(finish)
      endif
      if (mod(i,dnout)==0) then
         call check(x,v,passed)
         call write_out(time,xall,np)
      endif

   enddo

   print*,'-------------------------------------------------------------------'
   print*,'Finished.'
   print*,'-------------------------------------------------------------------'

end program test
