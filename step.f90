module step
 implicit none

contains
 subroutine step_leapfrog(ndim,x,v,gcov,gcon,fterm,dt)
  use metric, only: get_sourceterms
  integer, intent(in) :: ndim
  real(8), dimension(1+ndim,1+ndim), intent(in) :: gcov, gcon
  real(8), dimension(ndim), intent(in) :: fterm
  real(8), dimension(ndim), intent(out) :: x,v
  real(8), intent(in) :: dt

  pmom = pmom + 0.5*dt*fterm ! Half step in position
  call cons2prim(p,v) ! Get v(phalf,x0)

  ! Initial first order prediction for position (xstar)
  xprev = x
  x = x + dt*v
  ! Converge to x
  do while (max(abs(xprev-x))<1.d-6)
     xprev = x
     call cons2prim(x,vstar) ! Get v(phalf,xstar)=vstar
     x = xprev + 0.5*dt*(vstar - v)
  enddo

  pmom_prev = pmom
  pmom = pmom + 0.5*dt*fterm !pmom_star
  ! Converge to p
  do while (max(abs(pmom_prev-pmom)))
    pmom_prev = pmom
    call cons2prim(pmom,v) ! Get vstar from pmom_star
    call get_sourceterms(ndim,x,v,gcov,gcon,fterm_star) ! Get fterm(pmom_star,x1)=fterm_star !!This will need to be get_fores
    pmom = pmom + 0.5*dt*(fterm_star - fterm)
 enddo

 call prim2cons(x,v,pmom) ! Return primitives x1,v1

 end subroutine step_leapfrog
end module step
