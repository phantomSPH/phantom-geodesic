module step
  implicit none

  integer, parameter :: &
          ileapfrog = 1, &
          irk2      = 2, &
          ieuler    = 3, &
          iheuns    = 4, &
          ilnro5    = 5

  !-- Default timestepping method
  integer :: steptype = ilnro5

  contains

  character(len=10) function stepname(i)
  integer, intent(in) :: i

  select case(i)
  case(ileapfrog)
     stepname = 'Leapfrog'
  case(irk2)
     stepname = 'RK2'
  case(ieuler)
     stepname = 'Euler'
  case(iheuns)
     stepname = "Heun's"
  case(ilnro5)
     stepname = 'L&R05'
  end select

  end function stepname

  subroutine timestep_all(xall,vall,np,energy,angmom,dt,gravity_between_particles,mall)
    use step_old,     only: timestep
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

    angmom = 0.
    energy = 0.

    ! if (gravity_between_particles) then

       ! loop over all the particles and determine the force term, fterm for each
       do i = 1,np
         x = xall(:,i)
         v = vall(:,i)
         call get_sourceterms(x,v,ftermone)
         call get_newtonian_force_new(np,xall,ftermone,mall,i)
         fterm(:,i) = ftermone(:)
       enddo


       ! Next we call the modified timestepping alogorithm
       !call step_heuns_all(xall,vall,fterm,dt,np,mall)
       !call step_landr05_all(xall,vall,fterm,dt,np,mall)
       ! call step_1(xall,vall,fterm,dt,np,mall)

       print*,"steptype",steptype
       select case(steptype)
       case(ileapfrog)
          call step_leapfrog(xall,vall,fterm,dt,np,mall)
       case(irk2)
          call step_rk2(xall,vall,fterm,dt,np,mall)
       case(ieuler)
          call step_1(xall,vall,fterm,dt,np,mall)
       case(iheuns)
          call step_heuns(xall,vall,fterm,dt,np,mall)
       case(ilnro5)
          call step_landr05(xall,vall,fterm,dt,np,mall)
       end select

       do j = 1, np

         call get_ev(xall(:,j),vall(:,j),energy_i,angmom_i)
         energy = energy + energy_i
         angmom = angmom + angmom_i

       enddo

       call get_newtonian_energy(np,xall,vall,energy,mall)


    ! else
    !    !$omp parallel default(none) &
    !    !$omp shared(i,np,xall,vall,dt,dnout,dtout) &
    !    !$omp private(j,x,v,passed,energy_i,angmom_i) &
    !    !$omp reduction(+:energy,angmom)
    !    !$omp do
    !    do j=1,np
    !       x = xall(:,j)
    !       v = vall(:,j)
    !       call timestep(dt,x,v)
    !
    !       xall(:,j) = x
    !       vall(:,j) = v
    !
    !       call get_ev(x,v,energy_i,angmom_i)
    !
    !       energy = energy + energy_i
    !       angmom = angmom + angmom_i
    !    enddo
    !    !$omp enddo
    !    !$omp end parallel
    ! endif
  end subroutine timestep_all

  !----------------------------------------------------------------
  !+
  !  Euler Method (1st order)
  !+
  !----------------------------------------------------------------
  subroutine step_1(x,v,fterm,dt,np,mall)
   use cons2prim, only: get_p_from_v, get_v_from_p
   integer, intent(in) :: np
   real, dimension(3,np), intent(in) :: fterm
   real, dimension(3,np), intent(inout) :: x,v
   real, dimension(np), intent(in) :: mall
   real :: pmom(3,np)
   real, intent(in) :: dt
   integer :: i

   do i = 1,np
      call get_p_from_v(pmom(:,i),v(:,i),x(:,i))
      pmom(:,i) = pmom(:,i) + dt*fterm(:,i)
      x(:,i)    = x(:,i)    + dt*v(:,i)
     call get_v_from_p(pmom(:,i),v(:,i),x(:,i))
  enddo

end subroutine step_1

  !----------------------------------------------------------------
  !+
  !  Heuns's Method (2nd Order)
  !+
  !----------------------------------------------------------------
  subroutine step_heuns(x,v,fterm,dt,np,mall)
   use force_gr, only: get_sourceterms,get_newtonian_force_new
   use cons2prim, only: get_p_from_v, get_v_from_p
   integer, intent(in) :: np
   real, dimension(3,np), intent(inout) :: fterm
   real, dimension(3,np), intent(inout) :: x,v
   real, dimension(np), intent(in) :: mall
   real, dimension(3,np) :: pmom, fterm_old, v_old, pmom_guess, x_guess, fterm_new
   real, intent(in) :: dt
   integer :: i


   do i=1,np
      call get_p_from_v(pmom(:,i),v(:,i),x(:,i))
      v_old(:,i)  = v(:,i)
      fterm_old(:,i)  = fterm(:,i)
      pmom_guess(:,i) = pmom(:,i) + dt*fterm(:,i)
      x_guess(:,i)    = x(:,i)    + dt*v(:,i)
      call get_v_from_p(pmom_guess(:,i),v(:,i),x_guess(:,i))
      call get_sourceterms(x_guess(:,i),v(:,i),fterm_new(:,i))
      call get_newtonian_force_new(np,x,fterm_new(:,i),mall,i)
   enddo

   do i=1,np
      pmom(:,i) = pmom(:,i) + 0.5*dt*(fterm_new(:,i) + fterm_old(:,i) )
      x(:,i)    = x(:,i)    + 0.5*dt*(v(:,i) + v_old(:,i))
      call get_v_from_p(pmom(:,i),v(:,i),x(:,i))
      call get_sourceterms(x(:,i),v(:,i),fterm(:,i))
      call get_newtonian_force_new(np,x,fterm(:,i),mall,i)
   end do


 end subroutine step_heuns

!----------------------------------------------------------------
!+
!  Modified leapfrog (2nd order) from Leimkuhler & Reich (2005)
!+
!----------------------------------------------------------------
subroutine step_landr05(x,v,fterm,dt,np,mall)
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
    pmom(:,i) = pmom(:,i) + 0.5*dt*(fterm(:,i))  !pmom_star


    ! Converge to p
    do while (.not. converged_pmom .and. iterations_pmom < max_iterations)
       iterations_pmom = iterations_pmom + 1
       pmom_prev(:,i) = pmom(:,i)
       call get_v_from_p(pmom(:,i),v(:,i),x(:,i))                ! Get vstar from pmom_star
       call get_sourceterms(x(:,i),v(:,i),fterm_star(:,i))       ! Get fterm(pmom_star,x1)=fterm_star !!This will need to be get_forces
       call get_newtonian_force_new(np,x,fterm_star(:,i),mall,i)
       ! fterm_star(:,i) = fterm(:,i)
       pmom(:,i) = pmom_prev(:,i) + 0.5*dt*(fterm_star(:,i) - fterm(:,i))

       if (maxval(abs(pmom_prev(:,i)-pmom(:,i)))<=ptol) converged_pmom = .true.
       fterm(:,i) = fterm_star(:,i)
    enddo


    if (.not. converged_pmom) print*, 'WARNING: implicit timestep did not & converge! pmom-pmom_prev =',&
    &                                    pmom(:,i)-pmom_prev(:,i)
enddo
 do i = 1, np
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
    call get_newtonian_force_new(np,x,fterm(:,i),mall,i)
 enddo

 do i = 1,np
    pmom(:,i) = pmom(:,i) + 0.5*dt*fterm(:,i)  ! Half step in position
    call get_v_from_p(pmom(:,i),v(:,i),x(:,i))
 enddo

end subroutine step_landr05
!----------------------------------------------------------------
!+
!  Modified leapfrog (2nd order)
!+
!----------------------------------------------------------------
subroutine step_leapfrog(x,v,fterm,dt,np,mall)
 use force_gr, only: get_sourceterms,get_newtonian_force_new,get_newtonian_force
 use cons2prim, only: get_p_from_v, get_v_from_p

 integer, intent(in) :: np
 real, dimension(3,np), intent(inout) :: fterm
 real, dimension(3,np), intent(inout) :: x,v
 real, dimension(np), intent(in) :: mall
 real, intent(in) :: dt
 real, dimension(3,np) :: pmom, vstar, fterm_star, xprev, pmom_prev
 real :: xtol, ptol
 logical :: converged_x, converged_pmom
 integer :: iterations_x, iterations_pmom, i
 integer, parameter :: max_iterations = 100

 xtol = 1.e-15
 ptol = 1.e-15

do i = 1, np

  call get_p_from_v(pmom(:,i),v(:,i),x(:,i)) ! primitive to conservative

  pmom(:,i) = pmom(:,i) + 0.5*dt*fterm(:,i)  ! Half step in position
  call get_v_from_p(pmom(:,i),v(:,i),x(:,i)) ! Get v(phalf,x0)

  ! Initial first order prediction for position (xstar)
  x(:,i) = x(:,i) + dt*v(:,i)

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


enddo

do i = 1, np
  call get_newtonian_force_new(np,x,fterm_star(:,i),mall,i)
  pmom(:,i) = pmom(:,i) + 0.5*dt*fterm_star(:,i)  !pmom_star
  !Converge to p
  do while (.not. converged_pmom .and. iterations_pmom < max_iterations)
     iterations_pmom = iterations_pmom + 1
     pmom_prev(:,i) = pmom(:,i)
     call get_v_from_p(pmom(:,i),v(:,i),x(:,i))                         ! Get vstar from pmom_star
     !call get_sourceterms(x(:,i),v(:,i),fterm_star(:,i))                ! Get fterm(pmom_star,x1)=fterm_star !!This will need to be get_forces
     pmom(:,i) = pmom_prev(:,i) + 0.5*dt*(fterm_star(:,i) - fterm(:,i))
     if (maxval(abs(pmom_prev(:,i)-pmom(:,i)))<=ptol) converged_pmom = .true.
     fterm(:,i) = fterm_star(:,i)
  enddo
  if (.not. converged_pmom) print*, 'WARNING: implicit timestep did not & converge! pmom-pmom_prev =',&
  &                                    pmom(:,i)-pmom_prev(:,i)

  call get_v_from_p(pmom(:,i),v(:,i),x(:,i))
enddo
end subroutine step_leapfrog

!----------------------------------------------------------------
!+
!  RK2 Method (2nd Order)
!+
!----------------------------------------------------------------
subroutine step_rk2(x,v,fterm,dt,np,mall)
 use force_gr, only: get_sourceterms,get_newtonian_force_new
 use cons2prim, only: get_p_from_v, get_v_from_p

 integer, intent(in) :: np
 real, dimension(3,np), intent(inout) :: fterm
 real, dimension(3,np), intent(inout) :: x,v
 real, dimension(np), intent(in) :: mall
 real, intent(in) :: dt
 real :: pmomstar(3,np),xstar(3,np),vstar(3,np),pmom(3,np)
 integer :: i

 do i = 1, np
    call get_p_from_v(pmom(:,i),v(:,i),x(:,i))
    ! ENTRY
    xstar(:,i)    = x(:,i)    + 0.5*dt*v(:,i)
    pmomstar(:,i) = pmom(:,i) + 0.5*dt*fterm(:,i)
    call get_v_from_p(pmomstar(:,i),vstar(:,i),xstar(:,i))
 enddo

 do i = 1,np
    call get_sourceterms(xstar(:,i),vstar(:,i),fterm(:,i))
    call get_newtonian_force_new(np,x,fterm(:,i),mall,i)
 enddo

 do i = 1,np
    x(:,i)        = x(:,i)    + dt*vstar(:,i)
    pmom(:,i)     = pmom(:,i) + dt*fterm(:,i)
 enddo

 do i = 1,np
    ! EXIT
    call get_v_from_p(pmom(:,i),v(:,i),x(:,i))
    call get_sourceterms(x(:,i),v(:,i),fterm(:,i))
    call get_newtonian_force_new(np,x,fterm(:,i),mall,i)
 enddo

end subroutine step_rk2

end module step
