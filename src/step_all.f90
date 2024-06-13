module step_all
  implicit none
contains

  subroutine timestep_all(xall,vall,np,energy,angmom,dt)
    use step,         only: timestep
    use energies,     only: get_ev
    real, allocatable, dimension(:,:), intent(inout) :: xall,vall
    integer, intent(in) :: np
    real, intent(in)    :: dt
    real, dimension(3) :: x,v
    real, intent(out) :: energy,angmom

    real :: energy_i, angmom_i
    integer :: j

    angmom =0.
    energy =0.

    !$omp parallel default(none) &
    !$omp shared(i,np,xall,vall,dt,dnout,dtout) &
    !$omp private(j,x,v,passed,energy_i,angmom_i) &
    !$omp reduction(+:energy,angmom)
    !$omp do
    do j=1,np
      print*,np,"np"
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

  end subroutine timestep_all

end module step_all
