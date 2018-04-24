!----------------------------------------------------------------
!+
!  Evolve the position of a point particle in a metric (e.g. Schwarzschild)
!+
!----------------------------------------------------------------
program test
 use init,         only: initialise
 use metric,       only: metric_type
 use metric_tools, only: coordinate_sys
 use force_gr,     only: get_sourceterms
 use step,         only: timestep, stepname, steptype, ilnro5
 use energies,     only: get_ev
 use utils_gr,     only: get_rderivs
 use output,       only: write_out, write_ev, write_xyz, write_vxyz
 use checks,       only: check,sanity_checks
 use utils,        only: timer
 use prompting,    only: prompt
 implicit none

 real :: dt, tmax, dtout

 real, allocatable, dimension(:,:) :: xall,vall
 real, dimension(3) :: x,v
 integer :: np
 real    :: time, energy_init, angmom_init, energy, angmom, energy_i, angmom_i
 integer :: nsteps,i,j,dnout,dnout_ev
 logical :: passed,write_pos_vel
 real    :: start,finish,tminus,frac_done,twall_elapsed,twallmax_approx
 integer :: percentage,prev_percent

 ! Defaults
 steptype = ilnro5
 dt       = 2390./1.e4  !10^4 steps per orbit (precessing orbit)
 tmax     = 2390.*4
 dnout_ev = 30
 dtout    = -1.
 write_pos_vel = .true.

 print*,'-------------------------------------------------------------------'
 print*,'GR-TEST'
 print*,'-------------------------------------------------------------------'
 print*,               'Metric type       = ',metric_type
 print*,               'Coord. sys. type  = ',coordinate_sys

 ! Set particles and perform checks
 call initialise(xall,vall,np,energy,angmom)

 call prompt(" Enter step type (1 = Leapfrog  |  2 = RK2  |  3 = Euler  |  4 = Heun's  |  5 = L&R05) ",steptype)
 call prompt(" Enter dt  ",dt,0.)
 call prompt(" Enter tmax",tmax,dt)
 call prompt(" Enter dtout (-ve don't write dumps)",dtout)
 call prompt(" Write to ev file every how many steps? ",dnout_ev,0)
 if (np>1) then
    write_pos_vel = .false.
    call prompt(" Write positions and velocities to one file?",write_pos_vel)
 endif

 nsteps  = int(tmax/dt)
 dnout   = int(dtout/dt)

 print*,''
 print*,               'Timestepping used = ',trim(stepname(steptype))
 write(*,'(a,f10.2)') ' dt                = ',dt
 write(*,'(a,f10.2)') ' tmax              = ',tmax
 write(*,'(a,i10)')   ' dnout_ev          = ',dnout_ev

 print*,'-------------------------------------------------------------------'
 print*,' Press ENTER to start...'
 print*,'-------------------------------------------------------------------'
 read*
 print*,'Go...!'

 time     = 0.
 energy_init = energy
 angmom_init = angmom

 if (dtout>0) call write_out(time,xall,np)
 if (dnout_ev>0) then
    call write_ev(time,energy-energy_init,angmom-angmom_init)
    if (write_pos_vel) then
       call write_xyz(time,xall,np)
       call write_vxyz(time,vall,np)
    endif
    call timer(start)
 endif

 prev_percent = 0
 do i=1,nsteps
    angmom =0.
    energy =0.

    time = time + dt
    !$omp parallel default(none) &
    !$omp shared(i,np,xall,vall,dt,dnout,dtout) &
    !$omp private(j,x,v,passed,energy_i,angmom_i) &
    !$omp reduction(+:energy,angmom)
    !$omp do
    do j=1,np
       x = xall(:,j)
       v = vall(:,j)
       call timestep(dt,x,v)
       if (mod(i,dnout)==0 .and. dtout>0) call check(x,v,passed)
       xall(:,j) = x
       vall(:,j) = v
       call get_ev(x,v,energy_i,angmom_i)
       energy = energy + energy_i
       angmom = angmom + angmom_i
    enddo
    !$omp enddo
    !$omp end parallel

    if (mod(i,dnout)==0 .and. dtout>0) then
      call write_out(time,xall,np)
    endif

    if (dnout_ev>0) then
    if (mod(i,dnout_ev)==0) then
       call write_ev(time,energy-energy_init,angmom-angmom_init)
       if (write_pos_vel) then
          call write_xyz(time,xall,np)
          call write_vxyz(time,vall,np)
       endif
       call timer(finish)
       twall_elapsed = finish-start
       frac_done     = time/tmax
       percentage    = nint(frac_done*100.)
       twallmax_approx = twall_elapsed/frac_done
       tminus        = twallmax_approx - twall_elapsed !(1.-frac_done)/frac_done*twall_elapsed
       if (percentage == prev_percent + 1) then
          if (tminus<100.) then
             write(*,'(i4.1,a,f10.2,a,f10.2)') percentage,'%   t =', time,'  estimated t-minus   (s):',tminus
          else
             write(*,'(i4.1,a,f10.2,a,f10.2)') percentage,'%   t =', time,'  estimated t-minus (min):',tminus/60.
          endif
       endif
       prev_percent = percentage
    endif
    endif

 enddo

 print*,'-------------------------------------------------------------------'
 print*,'Finished.'
 print*,'-------------------------------------------------------------------'

end program test
