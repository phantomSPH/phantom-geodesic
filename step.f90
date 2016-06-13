module step
 implicit none

contains
 !----------------------------------------------------------------
 !+
 !  Modified leapfrom (2nd order)
 !+
 !----------------------------------------------------------------
 subroutine step_leapfrog(x,v,fterm,dt)
  use metric, only: get_sourceterms
  use cons2prim, only: get_p_from_v, get_v_from_p
  real, dimension(3), intent(in) :: fterm
  real, dimension(3), intent(out) :: x,v
  real, dimension(3) :: pmom, vstar, fterm_star, xprev, pmom_prev
  real, intent(in) :: dt
  real :: xtol, ptol
  integer :: iterations
  integer, parameter :: max_iterations = 5000

  iterations = 0
  xtol = 1.d-3
  ptol = 1.d-3

  call get_p_from_v(pmom,v,x) ! primitive to conservative
  pmom = pmom + 0.5*dt*fterm ! Half step in position
  call get_v_from_p(pmom,v,x) ! Get v(phalf,x0)
  ! Initial first order prediction for position (xstar)
  xprev = x
  x = x + dt*v
  ! Converge to x
  do while (maxval(abs(xprev-x))>=xtol .and. iterations < max_iterations)
   iterations = iterations + 1
   xprev = x
   call get_v_from_p(pmom,vstar,x) ! Get v(phalf,xstar)=vstar
   x = xprev + 0.5*dt*(vstar - v)
  enddo
  if (iterations >= max_iterations) print*, 'WARNING: implicit timestep did not &
  &converge! maxval(abs(xprev-x)) =', maxval(abs(xprev-x)), iterations
  iterations = 0
  pmom_prev = pmom
  pmom = pmom + 0.5*dt*fterm !pmom_star
  ! Converge to p
  do while (maxval(abs(pmom_prev-pmom))>=ptol .and. iterations < max_iterations)
   iterations = iterations + 1
   pmom_prev = pmom
   call get_v_from_p(pmom,v,x)               ! Get vstar from pmom_star
   call get_sourceterms(x,v,fterm_star) ! Get fterm(pmom_star,x1)=fterm_star !!This will need to be get_fores
   pmom = pmom_prev + 0.5*dt*(fterm_star - fterm)
  enddo
  if (iterations >= max_iterations) print*, 'WARNING: implicit timestep did not &
  &converge! pmom-pmom_prev =',pmom-pmom_prev
  call get_v_from_p(pmom,v,x)
 end subroutine step_leapfrog


 !----------------------------------------------------------------
 !+
 !  Euler Method (1st order)
 !+
 !----------------------------------------------------------------
 subroutine step_1(x,v,fterm,dt)
  use cons2prim, only: get_p_from_v, get_v_from_p
  real, dimension(3), intent(in) :: fterm
  real, dimension(3), intent(out) :: x,v
  real :: pmom(3)
  real, intent(in) :: dt

  call get_p_from_v(pmom,v,x)
  pmom = pmom + dt*fterm
  x    = x    + dt*v
  call get_v_from_p(pmom,v,x)
 end subroutine


 !----------------------------------------------------------------
 !+
 !  Heuns's Method (2nd Order)
 !+
 !----------------------------------------------------------------
 subroutine step_heuns(x,v,fterm,dt)
  use metric, only: get_sourceterms
  use cons2prim, only: get_p_from_v, get_v_from_p
  real, dimension(3), intent(in) :: fterm
  real, dimension(3), intent(out) :: x,v
  real, dimension(3) :: pmom, fterm_old, v_old, pmom_guess, x_guess, fterm_new
  real, intent(in) :: dt

  call get_p_from_v(pmom,v,x)
  v_old      = v
  fterm_old  = fterm
  pmom_guess = pmom + dt*fterm
  x_guess    = x    + dt*v
  call get_v_from_p(pmom_guess,v,x_guess)
  call get_sourceterms(x_guess,v,fterm_new)
  pmom = pmom + 0.5*dt*(fterm_new + fterm_old )
  x    = x    + 0.5*dt*(v + v_old)
  call get_v_from_p(pmom,v,x)

 end subroutine
end module step
