module options
implicit none
! 2390./1.e4
! 2390.*4
 !-- Defaults
 real    :: dt            = 10 !10^4 steps per orbit (precessing orbit)
 real    :: tmax          = 1e8
 integer :: dnout_ev      = 30
 real    :: dtout         = -1.
 logical :: write_pos_vel = .true.

end module options
