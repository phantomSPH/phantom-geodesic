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
 use step,         only: timestep, stepname, steptype, ilnro5
 use utils_gr,     only: get_ev
 use output,       only: write_out, write_ev, write_xyz, write_vxyz
 use checks,       only: check,sanity_checks
 use utils_gr,     only: get_rderivs
 use prompting,    only: prompt
 implicit none

 integer, parameter :: ndumps=15000
 real :: dt, tmax, dtout_ev, dtout

 real, allocatable, dimension(:,:) :: xall,vall
 real, dimension(3) :: x,v
 integer :: np
 real    :: time, energy_init, angmom_init, energy, angmom, energy_i, angmom_i
 integer :: nsteps, i,j, dnout, dnout_ev
 logical :: passed
 real    :: start, finish, frac_done, twall_elapsed
 real    :: tminus
 integer :: percentage,prev_percent

 steptype = ilnro5
 dt       = 1.e-2
 tmax     = 30000
 call prompt(" Enter step type (1 = Leapfrog  |  2 = RK2  |  3 = Euler  |  4 = Heun's  |  5 = L&R05) ",steptype)
 call prompt(" Enter dt  ",dt,0.)
 call prompt(" Enter tmax",tmax,dt)
 dtout_ev = tmax/ndumps
 dtout    = dtout_ev*1000.

 print*,'-------------------------------------------------------------------'
 print*,'GR-TEST'
 print*,'-------------------------------------------------------------------'
 print*,               'Metric type       = ',metric_type
 print*,               'Coord. sys. type  = ',coordinate_sys
 print*,               'Timestepping used = ',trim(stepname(steptype))
 write(*,'(a,f10.2)') ' dt                = ',dt
 write(*,'(a,f10.2)') ' tmax              = ',tmax
 write(*,'(a,f10.2)') ' dtout_ev          = ',dtout_ev
 nsteps   = int(tmax/dt)
 dnout    = int(dtout/dt)
 dnout_ev = int(dtout_ev/dt)
 time     = 0.
 call setup(xall, vall,np)
 print*,'-------------------------------------------------------------------'
 print*,' Press ENTER to start...'
 print*,'-------------------------------------------------------------------'
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
