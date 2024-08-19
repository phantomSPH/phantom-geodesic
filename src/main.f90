!----------------------------------------------------------------
!+
!  Evolve the position of a point particle in a metric (e.g. Schwarzschild)
!+
!----------------------------------------------------------------
program test
 use init,         only: initialise
 use metric,       only: metric_type,a
 use metric_tools, only: coordinate_sys
 use force_gr,     only: get_sourceterms
 use step,         only: stepname, steptype
 use utils_gr,     only: get_rderivs
 use output,       only: write_out, write_ev, write_xyz, write_vxyz
 use checks,       only: check,sanity_checks
 use utils,        only: timer
 use prompting,    only: prompt
 use options,      only: dt, tmax, dtout, dnout_ev, write_pos_vel
 use infile,       only: init_infile
 use step_all,     only: timestep_all
 implicit none

 real, allocatable, dimension(:,:) :: xall,vall
 real, allocatable, dimension(:)   :: mall
 integer :: np
 real    :: time, energy_init, angmom_init, energy, angmom
 integer :: nsteps,i,j,dnout
 logical :: passed
 real    :: start,finish,tminus,frac_done,twall_elapsed,twallmax_approx
 integer :: percentage,prev_percent
 logical :: use_self_gravity

 use_self_gravity = .true.

 print*,'-------------------------------------------------------------------'
 print*,'GR-TEST'
 print*,'-------------------------------------------------------------------'

 time = 0.

 ! Read options from infile
 call init_infile('grtest.in')

 print*,              'Metric type       = ',trim(metric_type)
 print*,              'Coord. sys. type  = ',trim(coordinate_sys)
 write(*,'(a,f5.3)') ' Black hole spin   = ',a

 ! Set particles and perform checks
 call initialise(time,xall,vall,np,energy,angmom,mall)
 nsteps  = int(tmax/dt)
 dnout   = int(dtout/dt)
 print*,''
 print*,               'Timestepping used = ',trim(stepname(steptype))
 write(*,'(a,f10.2)') ' dt                = ',dt
 write(*,'(a,f10.2)') ' tmax              = ',tmax
 write(*,'(a,i10)')   ' dnout_ev          = ',dnout_ev
 print*,'-------------------------------------------------------------------'
 print*,'Start:'
 print*,'-------------------------------------------------------------------'
 print*,              'Metric type       = ',trim(metric_type)
 energy_init = energy
 angmom_init = angmom
 if (dtout>0.) call write_out(time,xall,vall,np,mall)
 if (dnout_ev>0) then
    call write_ev(time,energy-energy_init,angmom-angmom_init,energy_init,angmom_init)
    if (write_pos_vel) then
       call write_xyz(time,xall,np)
       call write_vxyz(time,vall,np)
    endif
    call timer(start)
 endif
 dt = 100
 prev_percent = 0
 do i=1,nsteps
    time = time + dt
    print*,"dt ,",dt
    call timestep_all(xall,vall,np,energy,angmom,dt,use_self_gravity,mall)
    do j=1,np
       if (dtout>0. .and. mod(i,dnout)==0) call check(xall(:,j),vall(:,j),passed)
    enddo

    if (dtout>0. .and. mod(i,dnout)==0) then
      call write_out(time,xall,vall,np,mall)
    endif

    if (dnout_ev>0) then
    if (mod(i,dnout_ev)==0) then
       call write_ev(time,energy-energy_init,angmom-angmom_init,energy_init,angmom_init)
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
