module energies
 implicit none

contains

subroutine get_ev(x,v,energy,angmom)
 use metric, only: metric_type, rs
 use metric_tools, only: coordinate_sys,get_metric
 use cons2prim, only:get_p_from_v
 use utils_gr, only:get_u0
 use utils,    only:cross_product
 real, intent(in), dimension(3) :: x,v
 real, intent(out) :: energy, angmom
 real, dimension(3) :: pmom,fourvel_space,angmomg3
 real :: v4(0:3), gcov(0:3,0:3), gcon(0:3,0:3), sqrtg
 real :: r, U0
 integer, save :: i = 0
 character(len=*), parameter :: force_type = 'GR'

 call get_u0(x,v,U0)
 call get_metric(x,gcov,gcon,sqrtg)
 v4(0) = 1.
 v4(1:3) = v(1:3)
 fourvel_space = U0*v
 if (metric_type=='Schwarzschild') then
    if (coordinate_sys=='Cartesian') then
       r      = sqrt(dot_product(x,x))
       energy = (1. - rs/r)*U0
       angmom = (x(1)*v(2)-x(2)*v(1))*U0
    else if (coordinate_sys=='Spherical') then
       energy = (1. - rs/x(1))*U0
       angmom = x(1)**2*v(3)*U0
    endif
 else if (metric_type == 'Minkowski' .and. force_type == 'Newtonian') then
    energy = 0.5*dot_product(v,v)
    angmom = x(1)*v(2)-x(2)*v(1)
 else if (metric_type == 'Kerr') then
    energy = -U0*v4(0)*dot_product(gcov(0,:),v4(:))
    if (coordinate_sys=='Spherical') then
       call get_p_from_v(pmom,v,x)
       angmom = pmom(3)
    else if (coordinate_sys=='Cartesian') then
       call cross_product(x,fourvel_space,angmomg3)
       angmom = angmomg3(3)
    endif
 else
    if (i==0) then
       i = i+1
       print*,'WARNING: Energy and angular momentum are not being calculated for this metric. They will just be set a huge.'
       print*,'Continue?'
       read*
       energy = huge(energy)
       angmom = huge(energy)
    endif
 endif

end subroutine get_ev

end module energies
