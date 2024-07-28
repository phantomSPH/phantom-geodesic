module step_all
  implicit none
contains

  subroutine timestep_all(xall,vall,np,energy,angmom,dt,gravity_between_particles,mall)
    use step,         only: timestep
    use energies,     only: get_ev,get_newtonian_energy
    use force_gr,     only: get_sourceterms,get_newtonian_force,get_newtonian_force_new

    real, allocatable, dimension(:,:), intent(inout) :: xall,vall
    real, allocatable, dimension(:), intent(inout) :: mall
    integer, intent(in) :: np
    real, intent(in)    :: dt
    real, dimension(3) :: x,v,ftermone
    real, intent(out) :: energy,angmom
    logical :: gravity_between_particles
    real, dimension(3,np) :: fterm

    real :: energy_i, angmom_i
    integer :: j,i

    angmom =0.
    energy =0.

    if (gravity_between_particles) then

       ! loop over all the particles and determine the force term, fterm for each
       do i = 1,np
         x = xall(:,i)
         v = vall(:,i)
         call get_sourceterms(x,v,ftermone)
         fterm(:,i) = ftermone(:)
       enddo

       call get_newtonian_force(np,xall,fterm,mall)
       ! Next we call the modified timestepping alogorithm
       !call step_heuns_all(xall,vall,fterm,dt,np,mall)
       call step_landr05_all(xall,vall,fterm,dt,np,mall)

       do j = 1, np

         call get_ev(xall(:,j),vall(:,j),energy_i,angmom_i)
         energy = energy + energy_i
         angmom = angmom + angmom_i

       enddo

       !call get_COM_angular_mom(xall,vall,mall,angmom,np)

       call get_newtonian_energy(np,xall,vall,energy,mall)

    else
       !$omp parallel default(none) &
       !$omp shared(i,np,xall,vall,dt,dnout,dtout) &
       !$omp private(j,x,v,passed,energy_i,angmom_i) &
       !$omp reduction(+:energy,angmom)
       !$omp do
       do j=1,np
          x = xall(:,j)
          v = vall(:,j)
          call timestep(dt,x,v)

          xall(:,j) = x
          vall(:,j) = v

          call get_ev(x,v,energy_i,angmom_i)

          energy = energy + energy_i
          angmom = angmom + angmom_i
       enddo
       !$omp enddo
       !$omp end parallel
    endif
  end subroutine timestep_all

  !----------------------------------------------------------------
  !+
  !  Heuns's Method (2nd Order)
  !+
  !----------------------------------------------------------------
  subroutine step_heuns_all(x,v,fterm,dt,np,mall)
   use force_gr, only: get_sourceterms,get_newtonian_force
   use cons2prim, only: get_p_from_v, get_v_from_p
   integer, intent(in) :: np
   real, dimension(3,np), intent(inout) :: fterm
   real, dimension(3,np), intent(inout) :: x,v
   real, dimension(np), intent(in) :: mall
   real, dimension(3,np) :: pmom, fterm_old, v_old, pmom_guess, x_guess, fterm_new
   real, intent(in) :: dt
   integer :: i
    real :: angmom,rvcross(3),mag_rv
    angmom = 0.


   do i=1,np
      call get_p_from_v(pmom(:,i),v(:,i),x(:,i))
      !print*,v(:,i),"initial velocity", pmom(:,i),"initial momentum"
      v_old(:,i)  = v(:,i)
      fterm_old(:,i)  = fterm(:,i)
      pmom_guess(:,i) = pmom(:,i) + dt*fterm(:,i)
      x_guess(:,i)    = x(:,i)    + dt*v(:,i)
      call get_v_from_p(pmom_guess(:,i),v(:,i),x_guess(:,i))
      call get_sourceterms(x_guess(:,i),v(:,i),fterm_new(:,i))
   enddo

   call get_newtonian_force(np,x_guess,fterm_new,mall)
   !print*,fterm_new,"fterm_new"

   do i=1,np
      pmom(:,i) = pmom(:,i) + 0.5*dt*(fterm_new(:,i) + fterm_old(:,i) )
      x(:,i)    = x(:,i)    + 0.5*dt*(v(:,i) + v_old(:,i))
      call get_v_from_p(pmom(:,i),v(:,i),x(:,i))
      !print*,v(:,i),"final velocity", pmom(:,i),"final momentum"
      call get_sourceterms(x(:,i),v(:,i),fterm(:,i))
   end do

   call get_newtonian_force(np,x,fterm,mall)


 end subroutine step_heuns_all


!----------------------------------------------------------------
!+
!  Modified leapfrog (2nd order) from Leimkuhler & Reich (2005)
!+
!----------------------------------------------------------------
subroutine step_landr05_all(x,v,fterm,dt,np,mall)
 use force_gr, only: get_sourceterms,get_newtonian_force_new,get_newtonian_force
 use cons2prim, only: get_p_from_v, get_v_from_p

 integer, intent(in) :: np
 real, dimension(3,np), intent(inout) :: fterm
 real, dimension(3,np), intent(inout) :: x,v
 real, dimension(np), intent(in) :: mall
 real, intent(in) :: dt

 real, dimension(3,np) :: pmom, vstar, fterm_star, xprev, pmom_prev
 real :: xtol, ptol, tol
 logical :: converged_x, converged_pmom
 integer :: iterations_x, iterations_pmom, i
 integer, parameter :: max_iterations = 10000
 real :: angmom,rvcross(3),mag_rv,x_rel(3,2),v_rel(3,2)

 converged_x = .false.
 converged_pmom = .false.
 iterations_x = 0
 iterations_pmom = 0
 angmom = 0

 tol  = 1.e-15
 xtol = tol
 ptol = tol
 print*,"USING LR05"
 do i = 1, np
    converged_pmom = .false.
    iterations_pmom = 0

    call get_p_from_v(pmom(:,i),v(:,i),x(:,i)) ! primitive to conservative
    pmom(:,i) = pmom(:,i) + 0.5*dt*fterm(:,i)  !pmom_star

    ! Converge to p
    do while (.not. converged_pmom .and. iterations_pmom < max_iterations)
       iterations_pmom = iterations_pmom + 1
       pmom_prev(:,i) = pmom(:,i)
       call get_v_from_p(pmom(:,i),v(:,i),x(:,i))                ! Get vstar from pmom_star
       call get_sourceterms(x(:,i),v(:,i),fterm_star(:,i))       ! Get fterm(pmom_star,x1)=fterm_star !!This will need to be get_forces
       call get_newtonian_force_new(np,x,fterm_star(:,i),mall,i)

       pmom(:,i) = pmom_prev(:,i) + 0.5*dt*(fterm_star(:,i) - fterm(:,i))
       if (maxval(abs(pmom_prev(:,i)-pmom(:,i)))<=ptol) converged_pmom = .true.
       fterm(:,i) = fterm_star(:,i)
    enddo

    if (.not. converged_pmom) print*, 'WARNING: implicit timestep did not & converge! pmom-pmom_prev =',&
    &                                    pmom(:,i)-pmom_prev(:,i)
    call get_v_from_p(pmom(:,i),v(:,i),x(:,i)) ! Get v(phalf,x0)
    x(:,i) = x(:,i) + dt*v(:,i)
 enddo

 do i = 1, np
    converged_x = .false.
    iterations_x = 0
    ! Converge to x
    do while ( .not. converged_x .and. iterations_x < max_iterations)
       iterations_x = iterations_x + 1
       xprev(:,i) = x(:,i)
       call get_v_from_p(pmom(:,i),vstar(:,i),x(:,i)) ! Get v(phalf,xstar)=vstar
       x(:,i) = xprev(:,i) + 0.5*dt*(vstar(:,i) - v(:,i))
       if (maxval(abs(xprev(:,i)-x(:,i)))<=xtol) converged_x = .true.
       v(:,i) = vstar(:,i)
    enddo
    if (.not. converged_x) print*, 'WARNING: implicit timestep did not converge! maxval(abs(xprev-x)) =',&
    &                                 maxval(abs(xprev(:,i)-x(:,i))), iterations_x

    call get_sourceterms(x(:,i),v(:,i),fterm(:,i))
 enddo

 call get_newtonian_force(np,x,fterm,mall)

 do i = 1,np
    pmom(:,i) = pmom(:,i) + 0.5*dt*fterm(:,i)  ! Half step in position
    call get_v_from_p(pmom(:,i),v(:,i),x(:,i))
 enddo

end subroutine step_landr05_all

!----------------------------------------------------------------
!+
!  COM frame of binary
!+
!----------------------------------------------------------------
subroutine get_COM_angular_mom(x,v,mall,angmom,np)
  use utils,    only:cross_product

  integer, intent(in) :: np
  real, dimension(3,np), intent(in) :: x,v
  real, dimension(np), intent(in) :: mall
  real, intent(out) :: angmom

  real :: x_com(3), v_com(3), xrel(3), vrel(3), angmomg3(3)
  real :: angmom_i
  integer :: i

  x_com = (x(:,1) + x(:,2))*0.5
  v_com = (v(:,1) + v(:,2))*0.5

  ! angmom = 0.
  angmom_i = 0.
  angmomg3 = 0.

  do i = 1, np
    xrel(:) = x(:,i) - x_com
    vrel(:) = v(:,i) - v_com
    call cross_product(xrel,vrel,angmomg3)
    angmom_i = sqrt(dot_product(angmomg3,angmomg3))
    print*,angmom_i,"angmom of np",i
    angmom = angmom + angmom_i
    !print*,angmom,"angmom in the COM frame"
  enddo
  print*,angmom,"angmom total"

end subroutine get_COM_angular_mom
end module step_all
