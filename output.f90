module output
 implicit none
contains
 subroutine write_out(time,x,v,energy_init,angmom_init)
  use utils_gr, only: dot_product_gr
  use metric, only: get_metric, mass1
  real, intent(in) :: time, x(1:3), v(1:3),energy_init,angmom_init
  real, dimension(0:3,0:3) :: gcov, gcon
  real :: sqrtg, v4(0:3), energy, U0, angmom

  v4(0)   = 1.
  v4(1:3) = v(1:3)
  call get_metric(x,gcov,gcon,sqrtg)
  U0 = 1./sqrt(-dot_product_gr(v4,v4,gcov))

  write(1,*) x
  write(2,*) v
  ! Schwarzschild energy and angmom
  energy = (1. - 2*mass1/sqrt(dot_product(x,x)))/sqrt(-dot_product_gr(v4,v4,gcov))
  angmom = (x(1)*v(2)-x(2)*v(1))*U0
  write(3,*) time, energy-energy_init, angmom-angmom_init 

 end subroutine write_out
end module
