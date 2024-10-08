module force_gr
 implicit none

contains

subroutine get_forcegr(x,v,dens,u,p,fterm)
 use metric,       only: mass1
 real,    intent(in)  :: x(3),v(3),dens,u,p
 real,    intent(out) :: fterm(3)
 real :: r,r2

 if (mass1 == 0.) then
   STOP 'Mass of BH is zero. Change it for this calculation'
 endif
 
 r2 = dot_product(x,x)
 r  = sqrt(r2)
 fterm = -mass1*x/(r2*r)

end subroutine get_forcegr

! Wrapper routine to call get_forcegr for a test particle
subroutine get_sourceterms(x,v,fterm)
 real, intent(in)  :: x(3),v(3)
 real, intent(out) :: fterm(3)
 real :: dens,u,p

 P = 0.
 u = 0.
 dens = 1. ! this value does not matter (will cancel in the momentum equation)
 call get_forcegr(x,v,dens,u,p,fterm)
end subroutine get_sourceterms


end module force_gr
