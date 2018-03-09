module eos
 implicit none

 real, parameter :: gam = 5./3.

contains

subroutine get_pressure(P,dens,u)
 real, intent(out), dimension(:) :: P
 real, intent(in), dimension(:)  :: dens, u

 P= (gam-1.0)*dens*u

end subroutine

subroutine get_u(u,P,dens)
 real, intent(in)  :: dens,P
 real, intent(out) :: u

 u = P/((gam-1.)*dens)

end subroutine

subroutine get_enthalpy(enth,dens,P)
 real, intent(in)  :: dens,P
 real, intent(out) :: enth

 ! Needed in dust case when dens = 0 causes P/dens = NaN and therefore enth = NaN
 if(abs(p)<tiny(p)) then
    enth=1.
 else
    enth = 1.+P/dens*(gam/(gam-1.))
 endif

end subroutine

end module
